// Copyright (c) the JPEG XL Project Authors. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdint.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "build/third_party/gflags/include/gflags/gflags.h"
#include "fetch_encoded.h"
#include "gflags/gflags.h"
#include "jxl/codestream_header.h"
#include "jxl/color_encoding.h"
#include "jxl/encode.h"
#include "jxl/thread_parallel_runner.h"
#include "jxl/types.h"
#include "lib/jxl/base/file_io.h"
#include "lib/jxl/base/padded_bytes.h"

DEFINE_bool(boolean_test_variable, false,
            "Testing that adding a variable works.");

namespace {

// RAII-wraps the C-API encoder.
class ManagedJxlEncoder {
 public:
  explicit ManagedJxlEncoder(size_t num_worker_threads)
      : encoder_(JxlEncoderCreate(NULL)),
        encoder_frame_settings_(JxlEncoderFrameSettingsCreate(encoder_, NULL)) {
    if (num_worker_threads > 1) {
      parallel_runner_ = JxlThreadParallelRunnerCreate(
          /*memory_manager=*/nullptr, num_worker_threads);
    }
  }
  ~ManagedJxlEncoder() {
    if (parallel_runner_ != nullptr) {
      JxlThreadParallelRunnerDestroy(parallel_runner_);
    }
    JxlEncoderDestroy(encoder_);
    if (compressed_buffer_) {
      free(compressed_buffer_);
    }
  }

  JxlEncoder* encoder_;
  JxlEncoderFrameSettings* encoder_frame_settings_;
  uint8_t* compressed_buffer_ = nullptr;
  size_t compressed_buffer_size_ = 0;
  size_t compressed_buffer_used_ = 0;
  void* parallel_runner_ = nullptr;  // TODO(tfish): fix type.
};

}  // namespace

int main(int argc, char** argv) {
  gflags::SetUsageMessage(
      std::string("JPEG XL-encodes an image.  Sample usage:\n") +
      std::string(argv[0]) +
      std::string(" <source_image_filename> <target_image_filename>"));
  uint32_t version = JxlEncoderVersion();
  std::stringstream version_string;
  version_string << version / 1000000 << "." << (version / 1000) % 1000 << "."
              << version % 1000 << std::endl;

  gflags::SetVersionString(version_string.str());
  if (argc == 1) {
    gflags::ShowUsageWithFlags(argv[0]);
    exit(1);
  }
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  return EXIT_SUCCESS;
}
