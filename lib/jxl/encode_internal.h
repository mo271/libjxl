/* Copyright (c) the JPEG XL Project Authors. All rights reserved.
 *
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#ifndef LIB_JXL_ENCODE_INTERNAL_H_
#define LIB_JXL_ENCODE_INTERNAL_H_

#include <jxl/encode.h>
#include <jxl/memory_manager.h>
#include <jxl/parallel_runner.h>
#include <jxl/types.h>
#include <sys/types.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

#include "lib/jxl/base/data_parallel.h"
#include "lib/jxl/base/padded_bytes.h"
#include "lib/jxl/base/status.h"
#include "lib/jxl/enc_aux_out.h"
#include "lib/jxl/enc_fast_lossless.h"
#include "lib/jxl/enc_frame.h"
#include "lib/jxl/memory_manager_internal.h"

namespace jxl {

/* Frame index box 'jxli' will start with Varint() for
NF: has type Varint(): number of frames listed in the index.
TNUM: has type u32: numerator of tick unit.
TDEN: has type u32: denominator of tick unit. Value 0 means the file is
ill-formed. per frame i listed: OFFi: has type Varint(): offset of start byte of
this frame compared to start byte of previous frame from this index in the JPEG
XL codestream. For the first frame, this is the offset from the first byte of
the JPEG XL codestream. Ti: has type Varint(): duration in ticks between the
start of this frame and the start of the next frame in the index. If this is the
last frame in the index, this is the duration in ticks between the start of this
frame and the end of the stream. A tick lasts TNUM / TDEN seconds. Fi: has type
Varint(): amount of frames the next frame in the index occurs after this frame.
If this is the last frame in the index, this is the amount of frames after this
frame in the remainder of the stream. Only frames that are presented by the
decoder are counted for this purpose, this excludes frames that are not intended
for display but for compositing with other frames, such as frames that aren't
the last frame with a duration of 0 ticks.

All the frames listed in jxli are keyframes and the first frame is
present in the list.
There shall be either zero or one Frame Index boxes in a JPEG XL file.
The offsets OFFi per frame are given as bytes in the codestream, not as
bytes in the file format using the box structure. This means if JPEG XL Partial
Codestream boxes are used, the offset is counted within the concatenated
codestream, bytes from box headers or non-codestream boxes are not counted.
*/

typedef struct JxlEncoderFrameIndexBoxEntryStruct {
  bool to_be_indexed;
  uint32_t duration;
  uint64_t OFFi;
} JxlEncoderFrameIndexBoxEntry;

typedef struct JxlEncoderFrameIndexBoxStruct {
  // We always need to record the first frame entry, so presence of the
  // first entry alone is not an indication if it was requested to be
  // stored.
  bool index_box_requested_through_api = false;

  int64_t NF() const { return entries.size(); }
  bool StoreFrameIndexBox() {
    for (auto e : entries) {
      if (e.to_be_indexed) {
        return true;
      }
    }
    return false;
  }
  int32_t TNUM = 1;
  int32_t TDEN = 1000;

  std::vector<JxlEncoderFrameIndexBoxEntry> entries;

  // That way we can ensure that every index box will have the first frame.
  // If the API user decides to mark it as an indexed frame, we call
  // the AddFrame again, this time with requested.
  void AddFrame(uint64_t OFFi, uint32_t duration, bool to_be_indexed) {
    // We call AddFrame to every frame.
    // Recording the first frame is required by the standard.
    // Knowing the last frame is required, since the last indexed frame
    // needs to know how many frames until the end.
    // To be able to tell how many frames there are between each index
    // entry we just record every frame here.
    if (entries.size() == 1) {
      if (OFFi == entries[0].OFFi) {
        // API use for the first frame, let's clear the already recorded first
        // frame.
        entries.clear();
      }
    }
    JxlEncoderFrameIndexBoxEntry e;
    e.to_be_indexed = to_be_indexed;
    e.OFFi = OFFi;
    e.duration = duration;
    entries.push_back(e);
  }
} JxlEncoderFrameIndexBox;

// The encoder options (such as quality, compression speed, ...) for a single
// frame, but not encoder-wide options such as box-related options.
typedef struct JxlEncoderFrameSettingsValuesStruct {
  // lossless is a separate setting from cparams because it is a combination
  // setting that overrides multiple settings inside of cparams.
  bool lossless;
  CompressParams cparams;
  JxlFrameHeader header;
  std::vector<JxlBlendInfo> extra_channel_blend_info;
  std::string frame_name;
  JxlBitDepth image_bit_depth;
  bool frame_index_box = false;
  jxl::AuxOut* aux_out = nullptr;
} JxlEncoderFrameSettingsValues;

typedef std::array<uint8_t, 4> BoxType;

// Utility function that makes a BoxType from a string literal. The string must
// have 4 characters, a 5th null termination character is optional.
constexpr BoxType MakeBoxType(const char* type) {
  return BoxType(
      {{static_cast<uint8_t>(type[0]), static_cast<uint8_t>(type[1]),
        static_cast<uint8_t>(type[2]), static_cast<uint8_t>(type[3])}});
}

constexpr unsigned char kContainerHeader[] = {
    0,   0,   0, 0xc, 'J',  'X', 'L', ' ', 0xd, 0xa, 0x87,
    0xa, 0,   0, 0,   0x14, 'f', 't', 'y', 'p', 'j', 'x',
    'l', ' ', 0, 0,   0,    0,   'j', 'x', 'l', ' '};

constexpr unsigned char kLevelBoxHeader[] = {0, 0, 0, 0x9, 'j', 'x', 'l', 'l'};

struct JxlEncoderQueuedFrame {
  JxlEncoderFrameSettingsValues option_values;
  ImageBundle frame;
  std::vector<uint8_t> ec_initialized;
};

struct JxlEncoderQueuedBox {
  BoxType type;
  std::vector<uint8_t> contents;
  bool compress_box;
};

using FJXLFrameUniquePtr =
    std::unique_ptr<JxlFastLosslessFrameState,
                    decltype(&JxlFastLosslessFreeFrameState)>;

// Either a frame, or a box, not both.
// Can also be a FJXL frame.
struct JxlEncoderQueuedInput {
  explicit JxlEncoderQueuedInput(const JxlMemoryManager& memory_manager)
      : frame(nullptr, jxl::MemoryManagerDeleteHelper(&memory_manager)),
        box(nullptr, jxl::MemoryManagerDeleteHelper(&memory_manager)) {}
  MemoryManagerUniquePtr<JxlEncoderQueuedFrame> frame;
  MemoryManagerUniquePtr<JxlEncoderQueuedBox> box;
  FJXLFrameUniquePtr fast_lossless_frame = {nullptr,
                                            JxlFastLosslessFreeFrameState};
};

static constexpr size_t kSmallBoxHeaderSize = 8;
static constexpr size_t kLargeBoxHeaderSize = 16;
static constexpr size_t kLargeBoxContentSizeThreshold =
    0x100000000ull - kSmallBoxHeaderSize;

size_t WriteBoxHeader(const jxl::BoxType& type, size_t size, bool unbounded,
                      bool force_large_box, uint8_t* output);

// Appends a JXL container box header with given type, size, and unbounded
// properties to output.
template <typename T>
void AppendBoxHeader(const jxl::BoxType& type, size_t size, bool unbounded,
                     T* output) {
  size_t current_size = output->size();
  output->resize(current_size + kLargeBoxHeaderSize);
  size_t header_size =
      WriteBoxHeader(type, size, unbounded, /*force_large_box=*/false,
                     output->data() + current_size);
  output->resize(current_size + header_size);
}

}  // namespace jxl

class JxlEncoderOutputProcessorWrapper {
  friend class Buffer;
  class Buffer {
    public:
      size_t size() const { return size_; };
      uint8_t* data() { return data_; }

      Buffer(uint8_t* buffer, size_t size, size_t bytes_used,
             JxlEncoderOutputProcessorWrapper* wrapper)
          : data_(buffer),
            size_(size),
            bytes_used_(bytes_used),
            wrapper_(wrapper) {}
      ~Buffer() {
        release();
      }

      Buffer(const Buffer&) = delete;
      Buffer(Buffer&& other) noexcept
          : Buffer(other.data_, other.size_, other.bytes_used_,
                   other.wrapper_) {
      other.data_ = nullptr;
      other.size_ = 0;
      }

      void advance(size_t count) {
      JXL_ASSERT(count <= size_);
      data_ += count;
      size_ -= count;
      bytes_used_ += count;
      }

      void release() {
        if (this->data_) {
        wrapper_->ReleaseBuffer(bytes_used_);
        }
        data_ = nullptr;
        size_ = 0;
      }

      Buffer& operator=(const Buffer&) = delete;
      Buffer& operator=(Buffer&& other) noexcept {
        data_ = other.data_;
        size_ = other.size_;
        wrapper_ = other.wrapper_;
        return *this;
      }
    private:
      uint8_t* data_;
      size_t size_;
      size_t bytes_used_;
      JxlEncoderOutputProcessorWrapper* wrapper_;
  };

 public:
  JxlEncoderOutputProcessorWrapper() = default;
  explicit JxlEncoderOutputProcessorWrapper(JxlEncoderOutputProcessor processor)
      : external_output_processor_(
            jxl::make_unique<JxlEncoderOutputProcessor>(processor)) {}

  bool HasAvailOut() const { return avail_out_ != nullptr; }

  // Caller can never overwrite a previously-written buffer. Asking for a buffer
  // with `min_size` such that `position + min_size` overlaps with a
  // previously-written buffer is invalid.
  jxl::StatusOr<Buffer> GetBuffer(size_t min_size, size_t requested_size = 0) {
      if (stop_requested_) return jxl::StatusCode::kNotEnoughBytes;
      if (requested_size == 0) requested_size = min_size;
      JXL_ASSERT(requested_size >= min_size);

      // If we support seeking, output_position_ == position_.
      // Otherwise, output_position_ <= position_.
      JXL_ASSERT(output_position_ <= position_);
      size_t additional_size = position_ - output_position_;

      if (external_output_processor_) {
        size_t size = requested_size + additional_size;
        uint8_t* user_buffer =
            static_cast<uint8_t*>(external_output_processor_->get_buffer(
                external_output_processor_->opaque, &size));
        if (size == 0 || user_buffer == nullptr) {
        stop_requested_ = true;
        return jxl::StatusCode::kNotEnoughBytes;
        }
        if (size < min_size + additional_size) {
        external_output_processor_->release_buffer(
            external_output_processor_->opaque, 0);
        } else {
        internal_buffers_.emplace(position_, InternalBuffer());
        return Buffer(user_buffer + additional_size, size - additional_size, 0,
                      this);
        }
      } else {
        if (min_size + additional_size < *avail_out_) {
        internal_buffers_.emplace(position_, InternalBuffer());
        return Buffer(*next_out_ + additional_size,
                      *avail_out_ - additional_size, 0, this);
        }
      }

      // Otherwise, we need to allocate our own buffer.
      auto it = internal_buffers_.emplace(position_, InternalBuffer()).first;
      InternalBuffer& buffer = it->second;
      size_t alloc_size = requested_size;
      it++;
      if (it != internal_buffers_.end()) {
        alloc_size = std::min(alloc_size, it->first - position_);
        JXL_ASSERT(alloc_size >= min_size);
      }
      buffer.owned_data.resize(alloc_size);
      return Buffer(buffer.owned_data.data(), alloc_size, 0, this);
  }

  void Seek(size_t pos) {
      if (external_output_processor_ && external_output_processor_->seek) {
        external_output_processor_->seek(external_output_processor_->opaque,
                                         pos);
        output_position_ = pos;
      }
      JXL_ASSERT(pos >= watermark_position_);
      position_ = pos;
  }

  void SetWatermark(size_t pos) {
      if (external_output_processor_) {
        external_output_processor_->set_watermark(
            external_output_processor_->opaque, pos);
      }
      watermark_position_ = pos;
      FlushOutput();
  }

  size_t CurrentPosition() const { return position_; }

  bool SetAvailOut(uint8_t** next_out, size_t* avail_out) {
      if (external_output_processor_) return false;
      avail_out_ = avail_out;
      next_out_ = next_out;
      FlushOutput();
      return true;
  }

  bool WasStopRequested() const { return stop_requested_; }
  bool HasOutputToWrite() const { return output_position_ < watermark_position_; }

 private:
  void ReleaseBuffer(size_t bytes_used) {
      auto it = internal_buffers_.find(position_);
      if (bytes_used == 0) {
        internal_buffers_.erase(it);
        return;
      }
      JXL_ASSERT(it != internal_buffers_.end());
      it->second.written_bytes = bytes_used;
      if (external_output_processor_ && external_output_processor_->seek &&
          !it->second.owned_data.empty()) {
        external_output_processor_->seek(external_output_processor_->opaque,
                                         position_);
        output_position_ = position_;
        while (output_position_ < position_ + bytes_used) {
        size_t num_to_write = position_ + bytes_used - output_position_;
        if (!AppendBufferToExternalProcessor(
                it->second.owned_data.data() + output_position_ - position_,
                num_to_write)) {
          return;
        }
        }
      }
      position_ += bytes_used;

      it++;
      if (it != internal_buffers_.end()) {
        JXL_ASSERT(it->first >= position_);
      }
  }

  // Tries to write all the bytes up to the watermark position.
  void FlushOutput() {
      while (output_position_ < watermark_position_ &&
             (avail_out_ == nullptr || *avail_out_ > 0)) {
        JXL_ASSERT(!internal_buffers_.empty());
        auto it = internal_buffers_.begin();
        // If this fails, we are trying to move the watermark past data that was
        // not written yet. This is a library programming error.
        JXL_ASSERT(output_position_ >= it->first);
        JXL_ASSERT(it->second.written_bytes != 0);
        size_t buffer_last_byte = it->first + it->second.written_bytes;
        if (!it->second.owned_data.empty()) {
        size_t start_in_buffer = output_position_ - it->first;
        // Guaranteed by the invariant on `internal_buffers_`.
        JXL_ASSERT(buffer_last_byte > output_position_);
        size_t num_to_write =
            std::min(buffer_last_byte, watermark_position_) - output_position_;
        if (avail_out_ != nullptr) {
          size_t n = std::min(num_to_write, *avail_out_);
          memcpy(*next_out_, it->second.owned_data.data() + start_in_buffer, n);
          *avail_out_ -= n;
          *next_out_ += n;
          output_position_ += n;
        } else {
          if (!AppendBufferToExternalProcessor(
                  it->second.owned_data.data() + start_in_buffer,
                  num_to_write)) {
            return;
          }
        }
        }
        if (buffer_last_byte == output_position_) {
        internal_buffers_.erase(it);
        }
      }
  }

  bool AppendBufferToExternalProcessor(void* data, size_t count) {
      JXL_ASSERT(external_output_processor_);
      size_t n = count;
      void* user_buffer = external_output_processor_->get_buffer(
          external_output_processor_->opaque, &n);
      if (user_buffer == nullptr || n == 0) {
        stop_requested_ = true;
        return false;
      }
      n = std::min(n, count);
      memcpy(user_buffer, data, n);
      external_output_processor_->release_buffer(
          external_output_processor_->opaque, n);
      output_position_ += n;
      return true;
  }

  struct InternalBuffer {
      // Bytes in the range `[consumed_bytes, written_bytes)` need to be flushed
      // out.
      size_t written_bytes = 0;
      // If data has been buffered, it is stored in `owned_data`.
      jxl::PaddedBytes owned_data;
  };

  // Invariant: `internal_buffers_` does not contain chunks that are entirely
  // below the output position.
  std::map<size_t, InternalBuffer> internal_buffers_;

  uint8_t** next_out_ = nullptr;
  size_t* avail_out_ = nullptr;
  size_t position_ = 0;
  size_t watermark_position_ = 0;
  // Either the position of the `external_output_processor_` or the position
  // `next_out_` points to.
  size_t output_position_ = 0;

  bool stop_requested_ = false;

  std::unique_ptr<JxlEncoderOutputProcessor> external_output_processor_;
};

// Internal use only struct, can only be initialized correctly by
// JxlEncoderCreate.
struct JxlEncoderStruct {
  JxlEncoderError error = JxlEncoderError::JXL_ENC_ERR_OK;
  JxlMemoryManager memory_manager;
  jxl::MemoryManagerUniquePtr<jxl::ThreadPool> thread_pool{
      nullptr, jxl::MemoryManagerDeleteHelper(&memory_manager)};
  JxlCmsInterface cms;
  std::vector<jxl::MemoryManagerUniquePtr<JxlEncoderFrameSettings>>
      encoder_options;

  size_t num_queued_frames;
  size_t num_queued_boxes;
  std::vector<jxl::JxlEncoderQueuedInput> input_queue;
  JxlEncoderOutputProcessorWrapper output_processor;

  // How many codestream bytes have been written, i.e.,
  // content of jxlc and jxlp boxes. Frame index box jxli
  // requires position indices to point to codestream bytes,
  // so we need to keep track of the total of flushed or queue
  // codestream bytes. These bytes may be in a single jxlc box
  // or across multiple jxlp boxes.
  size_t codestream_bytes_written_beginning_of_frame;
  size_t codestream_bytes_written_end_of_frame;
  jxl::JxlEncoderFrameIndexBox frame_index_box;

  // Force using the container even if not needed
  bool use_container;
  // User declared they will add metadata boxes
  bool use_boxes;

  // TODO(lode): move level into jxl::CompressParams since some C++
  // implementation decisions should be based on it: level 10 allows more
  // features to be used.
  int32_t codestream_level;
  bool store_jpeg_metadata;
  jxl::CodecMetadata metadata;
  std::vector<uint8_t> jpeg_metadata;

  // Wrote any output at all, so wrote the data before the first user added
  // frame or box, such as signature, basic info, ICC profile or jpeg
  // reconstruction box.
  bool wrote_bytes;
  jxl::CompressParams last_used_cparams;
  JxlBasicInfo basic_info;

  // Encoder wrote a jxlp (partial codestream) box, so any next codestream
  // parts must also be written in jxlp boxes, a single jxlc box cannot be
  // used. The counter is used for the 4-byte jxlp box index header.
  size_t jxlp_counter;

  bool frames_closed;
  bool boxes_closed;
  bool basic_info_set;
  bool color_encoding_set;
  bool intensity_target_set;
  bool allow_expert_options = false;
  int brotli_effort = -1;

  // Takes the first frame in the input_queue, encodes it, and appends
  // the bytes to the output_byte_queue.
  jxl::Status RefillOutputByteQueue();

  bool MustUseContainer() const {
    return use_container || (codestream_level != 5 && codestream_level != -1) ||
           store_jpeg_metadata || use_boxes;
  }

  // `write_box` must never seek before the position the output wrapper was at
  // the moment of the call, and must leave the output wrapper such that its
  // position is one byte past the end of the written box.
  // TODO(veluca): add a function for appending a box that we already have the contents of.
  template <typename WriteBox>
  jxl::Status AppendBox(const jxl::BoxType& type, bool unbounded,
                        size_t box_max_size, const WriteBox& write_box);
};

struct JxlEncoderFrameSettingsStruct {
  JxlEncoder* enc;
  jxl::JxlEncoderFrameSettingsValues values;
};

struct JxlEncoderStatsStruct {
  jxl::AuxOut aux_out;
};

#endif  // LIB_JXL_ENCODE_INTERNAL_H_
