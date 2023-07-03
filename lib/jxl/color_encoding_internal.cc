// Copyright (c) the JPEG XL Project Authors. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "lib/jxl/color_encoding_internal.h"

#include <errno.h>

#include <array>
#include <cmath>

#include "lib/jxl/color_management.h"
#include "lib/jxl/common.h"
#include "lib/jxl/fields.h"
#include "lib/jxl/jxl_cms_internal.h"
#include "lib/jxl/matrix_ops.h"

#if JPEGXL_ENABLE_SKCMS
#include "lib/jxl/jxl_skcms.h"
#else  // JPEGXL_ENABLE_SKCMS
#include "lcms2.h"
#include "lcms2_plugin.h"
#endif  // JPEGXL_ENABLE_SKCMS


namespace jxl {
namespace {

// Highest reasonable value for the gamma of a transfer curve.
constexpr uint32_t kMaxGamma = 8192;

// These strings are baked into Description - do not change.

std::string ToString(ColorSpace color_space) {
  switch (color_space) {
    case ColorSpace::kRGB:
      return "RGB";
    case ColorSpace::kGray:
      return "Gra";
    case ColorSpace::kXYB:
      return "XYB";
    case ColorSpace::kUnknown:
      return "CS?";
  }
  // Should not happen - visitor fails if enum is invalid.
  JXL_ABORT("Invalid ColorSpace %u", static_cast<uint32_t>(color_space));
}

std::string ToString(WhitePoint white_point) {
  switch (white_point) {
    case WhitePoint::kD65:
      return "D65";
    case WhitePoint::kCustom:
      return "Cst";
    case WhitePoint::kE:
      return "EER";
    case WhitePoint::kDCI:
      return "DCI";
  }
  // Should not happen - visitor fails if enum is invalid.
  JXL_ABORT("Invalid WhitePoint %u", static_cast<uint32_t>(white_point));
}

std::string ToString(Primaries primaries) {
  switch (primaries) {
    case Primaries::kSRGB:
      return "SRG";
    case Primaries::k2100:
      return "202";
    case Primaries::kP3:
      return "DCI";
    case Primaries::kCustom:
      return "Cst";
  }
  // Should not happen - visitor fails if enum is invalid.
  JXL_ABORT("Invalid Primaries %u", static_cast<uint32_t>(primaries));
}

std::string ToString(TransferFunction transfer_function) {
  switch (transfer_function) {
    case TransferFunction::kSRGB:
      return "SRG";
    case TransferFunction::kLinear:
      return "Lin";
    case TransferFunction::k709:
      return "709";
    case TransferFunction::kPQ:
      return "PeQ";
    case TransferFunction::kHLG:
      return "HLG";
    case TransferFunction::kDCI:
      return "DCI";
    case TransferFunction::kUnknown:
      return "TF?";
  }
  // Should not happen - visitor fails if enum is invalid.
  JXL_ABORT("Invalid TransferFunction %u",
            static_cast<uint32_t>(transfer_function));
}

std::string ToString(RenderingIntent rendering_intent) {
  switch (rendering_intent) {
    case RenderingIntent::kPerceptual:
      return "Per";
    case RenderingIntent::kRelative:
      return "Rel";
    case RenderingIntent::kSaturation:
      return "Sat";
    case RenderingIntent::kAbsolute:
      return "Abs";
  }
  // Should not happen - visitor fails if enum is invalid.
  JXL_ABORT("Invalid RenderingIntent %u",
            static_cast<uint32_t>(rendering_intent));
}

static double F64FromCustomxyI32(const int32_t i) { return i * 1E-6; }
static Status F64ToCustomxyI32(const double f, int32_t* JXL_RESTRICT i) {
  if (!(-4 <= f && f <= 4)) {
    return JXL_FAILURE("F64 out of bounds for CustomxyI32");
  }
  *i = static_cast<int32_t>(roundf(f * 1E6));
  return true;
}

Status ConvertExternalToInternalWhitePoint(const JxlWhitePoint external,
                                           WhitePoint* internal) {
  switch (external) {
    case JXL_WHITE_POINT_D65:
      *internal = WhitePoint::kD65;
      return true;
    case JXL_WHITE_POINT_CUSTOM:
      *internal = WhitePoint::kCustom;
      return true;
    case JXL_WHITE_POINT_E:
      *internal = WhitePoint::kE;
      return true;
    case JXL_WHITE_POINT_DCI:
      *internal = WhitePoint::kDCI;
      return true;
  }
  return JXL_FAILURE("Invalid WhitePoint enum value");
}

Status ConvertExternalToInternalPrimaries(const JxlPrimaries external,
                                          Primaries* internal) {
  switch (external) {
    case JXL_PRIMARIES_SRGB:
      *internal = Primaries::kSRGB;
      return true;
    case JXL_PRIMARIES_CUSTOM:
      *internal = Primaries::kCustom;
      return true;
    case JXL_PRIMARIES_2100:
      *internal = Primaries::k2100;
      return true;
    case JXL_PRIMARIES_P3:
      *internal = Primaries::kP3;
      return true;
  }
  return JXL_FAILURE("Invalid Primaries enum value");
}

Status ConvertExternalToInternalTransferFunction(
    const JxlTransferFunction external, TransferFunction* internal) {
  switch (external) {
    case JXL_TRANSFER_FUNCTION_709:
      *internal = TransferFunction::k709;
      return true;
    case JXL_TRANSFER_FUNCTION_UNKNOWN:
      *internal = TransferFunction::kUnknown;
      return true;
    case JXL_TRANSFER_FUNCTION_LINEAR:
      *internal = TransferFunction::kLinear;
      return true;
    case JXL_TRANSFER_FUNCTION_SRGB:
      *internal = TransferFunction::kSRGB;
      return true;
    case JXL_TRANSFER_FUNCTION_PQ:
      *internal = TransferFunction::kPQ;
      return true;
    case JXL_TRANSFER_FUNCTION_DCI:
      *internal = TransferFunction::kDCI;
      return true;
    case JXL_TRANSFER_FUNCTION_HLG:
      *internal = TransferFunction::kHLG;
      return true;
    case JXL_TRANSFER_FUNCTION_GAMMA:
      return JXL_FAILURE("Gamma should be handled separately");
  }
  return JXL_FAILURE("Invalid TransferFunction enum value");
}

Status ConvertExternalToInternalRenderingIntent(
    const JxlRenderingIntent external, RenderingIntent* internal) {
  switch (external) {
    case JXL_RENDERING_INTENT_PERCEPTUAL:
      *internal = RenderingIntent::kPerceptual;
      return true;
    case JXL_RENDERING_INTENT_RELATIVE:
      *internal = RenderingIntent::kRelative;
      return true;
    case JXL_RENDERING_INTENT_SATURATION:
      *internal = RenderingIntent::kSaturation;
      return true;
    case JXL_RENDERING_INTENT_ABSOLUTE:
      *internal = RenderingIntent::kAbsolute;
      return true;
  }
  return JXL_FAILURE("Invalid RenderingIntent enum value");
}

bool ApplyCICP(const uint8_t color_primaries,
               const uint8_t transfer_characteristics,
               const uint8_t matrix_coefficients, const uint8_t full_range,
               ColorEncoding* JXL_RESTRICT c) {
  if (matrix_coefficients != 0) return false;
  if (full_range != 1) return false;

  const auto primaries = static_cast<Primaries>(color_primaries);
  const auto tf = static_cast<TransferFunction>(transfer_characteristics);
  if (tf == TransferFunction::kUnknown || !EnumValid(tf)) return false;
  if (primaries == Primaries::kCustom ||
      !(color_primaries == 12 || EnumValid(primaries))) {
    return false;
  }
  c->SetColorSpace(ColorSpace::kRGB);
  c->tf.SetTransferFunction(tf);
  if (primaries == Primaries::kP3) {
    c->white_point = WhitePoint::kDCI;
    c->primaries = Primaries::kP3;
  } else if (color_primaries == 12) {
    c->white_point = WhitePoint::kD65;
    c->primaries = Primaries::kP3;
  } else {
    c->white_point = WhitePoint::kD65;
    c->primaries = primaries;
  }
  return true;
}


#if JPEGXL_ENABLE_SKCMS

ColorSpace ColorSpaceFromProfile(const skcms_ICCProfile& profile) {
  switch (profile.data_color_space) {
    case skcms_Signature_RGB:
    case skcms_Signature_CMYK:
      // spec says CMYK is encoded as RGB (the kBlack extra channel signals that
      // it is actually CMYK)
      return ColorSpace::kRGB;
    case skcms_Signature_Gray:
      return ColorSpace::kGray;
    default:
      return ColorSpace::kUnknown;
  }
}

// vector_out := matmul(matrix, vector_in)
void MatrixProduct(const skcms_Matrix3x3& matrix, const float vector_in[3],
                   float vector_out[3]) {
  for (int i = 0; i < 3; ++i) {
    vector_out[i] = 0;
    for (int j = 0; j < 3; ++j) {
      vector_out[i] += matrix.vals[i][j] * vector_in[j];
    }
  }
}


#else  // JPEGXL_ENABLE_SKCMS

uint32_t Type32(const ColorEncoding& c, bool cmyk) {
  if (cmyk) return TYPE_CMYK_FLT;
  if (c.IsGray()) return TYPE_GRAY_FLT;
  return TYPE_RGB_FLT;
}

uint32_t Type64(const ColorEncoding& c) {
  if (c.IsGray()) return TYPE_GRAY_DBL;
  return TYPE_RGB_DBL;
}

ColorSpace ColorSpaceFromProfile(const Profile& profile) {
  switch (cmsGetColorSpace(profile.get())) {
    case cmsSigRgbData:
    case cmsSigCmykData:
      return ColorSpace::kRGB;
    case cmsSigGrayData:
      return ColorSpace::kGray;
    default:
      return ColorSpace::kUnknown;
  }
}

// "profile1" is pre-decoded to save time in DetectTransferFunction.
Status ProfileEquivalentToICC(const cmsContext context, const Profile& profile1,
                              const PaddedBytes& icc, const ColorEncoding& c) {
  const uint32_t type_src = Type64(c);

  Profile profile2;
  JXL_RETURN_IF_ERROR(DecodeProfile(context, icc, &profile2));

  Profile profile_xyz;
  JXL_RETURN_IF_ERROR(CreateProfileXYZ(context, &profile_xyz));

  const uint32_t intent = INTENT_RELATIVE_COLORIMETRIC;
  const uint32_t flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_BLACKPOINTCOMPENSATION |
                         cmsFLAGS_HIGHRESPRECALC;
  Transform xform1(cmsCreateTransformTHR(context, profile1.get(), type_src,
                                         profile_xyz.get(), TYPE_XYZ_DBL,
                                         intent, flags));
  Transform xform2(cmsCreateTransformTHR(context, profile2.get(), type_src,
                                         profile_xyz.get(), TYPE_XYZ_DBL,
                                         intent, flags));
  if (xform1 == nullptr || xform2 == nullptr) {
    return JXL_FAILURE("Failed to create transform");
  }

  double in[3];
  double out1[3];
  double out2[3];

  // Uniformly spaced samples from very dark to almost fully bright.
  const double init = 1E-3;
  const double step = 0.2;

  if (c.IsGray()) {
    // Finer sampling and replicate each component.
    for (in[0] = init; in[0] < 1.0; in[0] += step / 8) {
      cmsDoTransform(xform1.get(), in, out1, 1);
      cmsDoTransform(xform2.get(), in, out2, 1);
      if (!ApproxEq(out1[0], out2[0], 2E-4)) {
        return false;
      }
    }
  } else {
    for (in[0] = init; in[0] < 1.0; in[0] += step) {
      for (in[1] = init; in[1] < 1.0; in[1] += step) {
        for (in[2] = init; in[2] < 1.0; in[2] += step) {
          cmsDoTransform(xform1.get(), in, out1, 1);
          cmsDoTransform(xform2.get(), in, out2, 1);
          for (size_t i = 0; i < 3; ++i) {
            if (!ApproxEq(out1[i], out2[i], 2E-4)) {
              return false;
            }
          }
        }
      }
    }
  }

  return true;
}

#endif  // JPEGXL_ENABLE_SKCMS

}  // namespace

CIExy Customxy::Get() const {
  CIExy xy;
  xy.x = F64FromCustomxyI32(x);
  xy.y = F64FromCustomxyI32(y);
  return xy;
}

Status Customxy::Set(const CIExy& xy) {
  JXL_RETURN_IF_ERROR(F64ToCustomxyI32(xy.x, &x));
  JXL_RETURN_IF_ERROR(F64ToCustomxyI32(xy.y, &y));
  size_t extension_bits, total_bits;
  if (!Bundle::CanEncode(*this, &extension_bits, &total_bits)) {
    return JXL_FAILURE("Unable to encode XY %f %f", xy.x, xy.y);
  }
  return true;
}

bool CustomTransferFunction::SetImplicit() {
  if (nonserialized_color_space == ColorSpace::kXYB) {
    if (!SetGamma(1.0 / 3)) JXL_ASSERT(false);
    return true;
  }
  return false;
}

Status CustomTransferFunction::SetGamma(double gamma) {
  if (gamma < (1.0f / kMaxGamma) || gamma > 1.0) {
    return JXL_FAILURE("Invalid gamma %f", gamma);
  }

  have_gamma_ = false;
  if (ApproxEq(gamma, 1.0)) {
    transfer_function_ = TransferFunction::kLinear;
    return true;
  }
  if (ApproxEq(gamma, 1.0 / 2.6)) {
    transfer_function_ = TransferFunction::kDCI;
    return true;
  }
  // Don't translate 0.45.. to kSRGB nor k709 - that might change pixel
  // values because those curves also have a linear part.

  have_gamma_ = true;
  gamma_ = roundf(gamma * kGammaMul);
  transfer_function_ = TransferFunction::kUnknown;
  return true;
}

namespace {

std::array<ColorEncoding, 2> CreateC2(const Primaries pr,
                                      const TransferFunction tf) {
  std::array<ColorEncoding, 2> c2;

  {
    ColorEncoding* c_rgb = c2.data() + 0;
    c_rgb->SetColorSpace(ColorSpace::kRGB);
    c_rgb->white_point = WhitePoint::kD65;
    c_rgb->primaries = pr;
    c_rgb->tf.SetTransferFunction(tf);
    JXL_CHECK(c_rgb->CreateICC());
  }

  {
    ColorEncoding* c_gray = c2.data() + 1;
    c_gray->SetColorSpace(ColorSpace::kGray);
    c_gray->white_point = WhitePoint::kD65;
    c_gray->primaries = pr;
    c_gray->tf.SetTransferFunction(tf);
    JXL_CHECK(c_gray->CreateICC());
  }

  return c2;
}

}  // namespace

const ColorEncoding& ColorEncoding::SRGB(bool is_gray) {
  static std::array<ColorEncoding, 2> c2 =
      CreateC2(Primaries::kSRGB, TransferFunction::kSRGB);
  return c2[is_gray];
}
const ColorEncoding& ColorEncoding::LinearSRGB(bool is_gray) {
  static std::array<ColorEncoding, 2> c2 =
      CreateC2(Primaries::kSRGB, TransferFunction::kLinear);
  return c2[is_gray];
}

CIExy ColorEncoding::GetWhitePoint() const {
  JXL_DASSERT(have_fields_);
  CIExy xy;
  switch (white_point) {
    case WhitePoint::kCustom:
      return white_.Get();

    case WhitePoint::kD65:
      xy.x = 0.3127;
      xy.y = 0.3290;
      return xy;

    case WhitePoint::kDCI:
      // From https://ieeexplore.ieee.org/document/7290729 C.2 page 11
      xy.x = 0.314;
      xy.y = 0.351;
      return xy;

    case WhitePoint::kE:
      xy.x = xy.y = 1.0 / 3;
      return xy;
  }
  JXL_ABORT("Invalid WhitePoint %u", static_cast<uint32_t>(white_point));
}

Status ColorEncoding::SetWhitePoint(const CIExy& xy) {
  JXL_DASSERT(have_fields_);
  if (xy.x == 0.0 || xy.y == 0.0) {
    return JXL_FAILURE("Invalid white point %f %f", xy.x, xy.y);
  }
  if (ApproxEq(xy.x, 0.3127) && ApproxEq(xy.y, 0.3290)) {
    white_point = WhitePoint::kD65;
    return true;
  }
  if (ApproxEq(xy.x, 1.0 / 3) && ApproxEq(xy.y, 1.0 / 3)) {
    white_point = WhitePoint::kE;
    return true;
  }
  if (ApproxEq(xy.x, 0.314) && ApproxEq(xy.y, 0.351)) {
    white_point = WhitePoint::kDCI;
    return true;
  }
  white_point = WhitePoint::kCustom;
  return white_.Set(xy);
}

PrimariesCIExy ColorEncoding::GetPrimaries() const {
  JXL_DASSERT(have_fields_);
  JXL_ASSERT(HasPrimaries());
  PrimariesCIExy xy;
  switch (primaries) {
    case Primaries::kCustom:
      xy.r = red_.Get();
      xy.g = green_.Get();
      xy.b = blue_.Get();
      return xy;

    case Primaries::kSRGB:
      xy.r.x = 0.639998686;
      xy.r.y = 0.330010138;
      xy.g.x = 0.300003784;
      xy.g.y = 0.600003357;
      xy.b.x = 0.150002046;
      xy.b.y = 0.059997204;
      return xy;

    case Primaries::k2100:
      xy.r.x = 0.708;
      xy.r.y = 0.292;
      xy.g.x = 0.170;
      xy.g.y = 0.797;
      xy.b.x = 0.131;
      xy.b.y = 0.046;
      return xy;

    case Primaries::kP3:
      xy.r.x = 0.680;
      xy.r.y = 0.320;
      xy.g.x = 0.265;
      xy.g.y = 0.690;
      xy.b.x = 0.150;
      xy.b.y = 0.060;
      return xy;
  }
  JXL_ABORT("Invalid Primaries %u", static_cast<uint32_t>(primaries));
}

Status ColorEncoding::SetPrimaries(const PrimariesCIExy& xy) {
  JXL_DASSERT(have_fields_);
  JXL_ASSERT(HasPrimaries());
  if (xy.r.x == 0.0 || xy.r.y == 0.0 || xy.g.x == 0.0 || xy.g.y == 0.0 ||
      xy.b.x == 0.0 || xy.b.y == 0.0) {
    return JXL_FAILURE("Invalid primaries %f %f %f %f %f %f", xy.r.x, xy.r.y,
                       xy.g.x, xy.g.y, xy.b.x, xy.b.y);
  }

  if (ApproxEq(xy.r.x, 0.64) && ApproxEq(xy.r.y, 0.33) &&
      ApproxEq(xy.g.x, 0.30) && ApproxEq(xy.g.y, 0.60) &&
      ApproxEq(xy.b.x, 0.15) && ApproxEq(xy.b.y, 0.06)) {
    primaries = Primaries::kSRGB;
    return true;
  }

  if (ApproxEq(xy.r.x, 0.708) && ApproxEq(xy.r.y, 0.292) &&
      ApproxEq(xy.g.x, 0.170) && ApproxEq(xy.g.y, 0.797) &&
      ApproxEq(xy.b.x, 0.131) && ApproxEq(xy.b.y, 0.046)) {
    primaries = Primaries::k2100;
    return true;
  }
  if (ApproxEq(xy.r.x, 0.680) && ApproxEq(xy.r.y, 0.320) &&
      ApproxEq(xy.g.x, 0.265) && ApproxEq(xy.g.y, 0.690) &&
      ApproxEq(xy.b.x, 0.150) && ApproxEq(xy.b.y, 0.060)) {
    primaries = Primaries::kP3;
    return true;
  }

  primaries = Primaries::kCustom;
  JXL_RETURN_IF_ERROR(red_.Set(xy.r));
  JXL_RETURN_IF_ERROR(green_.Set(xy.g));
  JXL_RETURN_IF_ERROR(blue_.Set(xy.b));
  return true;
}

Status ColorEncoding::CreateICC() {
  InternalRemoveICC();
  if (!MaybeCreateProfile(*this, &icc_)) {
    return JXL_FAILURE("Failed to create profile from fields");
  }
  return true;
}

std::string Description(const ColorEncoding& c_in) {
  // Copy required for Implicit*
  ColorEncoding c = c_in;

  std::string d = ToString(c.GetColorSpace());

  if (!c.ImplicitWhitePoint()) {
    d += '_';
    if (c.white_point == WhitePoint::kCustom) {
      const CIExy wp = c.GetWhitePoint();
      d += ToString(wp.x) + ';';
      d += ToString(wp.y);
    } else {
      d += ToString(c.white_point);
    }
  }

  if (c.HasPrimaries()) {
    d += '_';
    if (c.primaries == Primaries::kCustom) {
      const PrimariesCIExy pr = c.GetPrimaries();
      d += ToString(pr.r.x) + ';';
      d += ToString(pr.r.y) + ';';
      d += ToString(pr.g.x) + ';';
      d += ToString(pr.g.y) + ';';
      d += ToString(pr.b.x) + ';';
      d += ToString(pr.b.y);
    } else {
      d += ToString(c.primaries);
    }
  }

  d += '_';
  d += ToString(c.rendering_intent);

  if (!c.tf.SetImplicit()) {
    d += '_';
    if (c.tf.IsGamma()) {
      d += 'g';
      d += ToString(c.tf.GetGamma());
    } else {
      d += ToString(c.tf.GetTransferFunction());
    }
  }

  return d;
}

Customxy::Customxy() { Bundle::Init(this); }
Status Customxy::VisitFields(Visitor* JXL_RESTRICT visitor) {
  uint32_t ux = PackSigned(x);
  JXL_QUIET_RETURN_IF_ERROR(visitor->U32(Bits(19), BitsOffset(19, 524288),
                                         BitsOffset(20, 1048576),
                                         BitsOffset(21, 2097152), 0, &ux));
  x = UnpackSigned(ux);
  uint32_t uy = PackSigned(y);
  JXL_QUIET_RETURN_IF_ERROR(visitor->U32(Bits(19), BitsOffset(19, 524288),
                                         BitsOffset(20, 1048576),
                                         BitsOffset(21, 2097152), 0, &uy));
  y = UnpackSigned(uy);
  return true;
}

CustomTransferFunction::CustomTransferFunction() { Bundle::Init(this); }
Status CustomTransferFunction::VisitFields(Visitor* JXL_RESTRICT visitor) {
  if (visitor->Conditional(!SetImplicit())) {
    JXL_QUIET_RETURN_IF_ERROR(visitor->Bool(false, &have_gamma_));

    if (visitor->Conditional(have_gamma_)) {
      // Gamma is represented as a 24-bit int, the exponent used is
      // gamma_ / 1e7. Valid values are (0, 1]. On the low end side, we also
      // limit it to kMaxGamma/1e7.
      JXL_QUIET_RETURN_IF_ERROR(visitor->Bits(24, kGammaMul, &gamma_));
      if (gamma_ > kGammaMul ||
          static_cast<uint64_t>(gamma_) * kMaxGamma < kGammaMul) {
        return JXL_FAILURE("Invalid gamma %u", gamma_);
      }
    }

    if (visitor->Conditional(!have_gamma_)) {
      JXL_QUIET_RETURN_IF_ERROR(
          visitor->Enum(TransferFunction::kSRGB, &transfer_function_));
    }
  }

  return true;
}

ColorEncoding::ColorEncoding() { Bundle::Init(this); }
Status ColorEncoding::VisitFields(Visitor* JXL_RESTRICT visitor) {
  if (visitor->AllDefault(*this, &all_default)) {
    // Overwrite all serialized fields, but not any nonserialized_*.
    visitor->SetDefault(this);
    return true;
  }

  JXL_QUIET_RETURN_IF_ERROR(visitor->Bool(false, &want_icc_));

  // Always send even if want_icc_ because this affects decoding.
  // We can skip the white point/primaries because they do not.
  JXL_QUIET_RETURN_IF_ERROR(visitor->Enum(ColorSpace::kRGB, &color_space_));

  if (visitor->Conditional(!WantICC())) {
    // Serialize enums. NOTE: we set the defaults to the most common values so
    // ImageMetadata.all_default is true in the common case.

    if (visitor->Conditional(!ImplicitWhitePoint())) {
      JXL_QUIET_RETURN_IF_ERROR(visitor->Enum(WhitePoint::kD65, &white_point));
      if (visitor->Conditional(white_point == WhitePoint::kCustom)) {
        JXL_QUIET_RETURN_IF_ERROR(visitor->VisitNested(&white_));
      }
    }

    if (visitor->Conditional(HasPrimaries())) {
      JXL_QUIET_RETURN_IF_ERROR(visitor->Enum(Primaries::kSRGB, &primaries));
      if (visitor->Conditional(primaries == Primaries::kCustom)) {
        JXL_QUIET_RETURN_IF_ERROR(visitor->VisitNested(&red_));
        JXL_QUIET_RETURN_IF_ERROR(visitor->VisitNested(&green_));
        JXL_QUIET_RETURN_IF_ERROR(visitor->VisitNested(&blue_));
      }
    }

    JXL_QUIET_RETURN_IF_ERROR(visitor->VisitNested(&tf));

    JXL_QUIET_RETURN_IF_ERROR(
        visitor->Enum(RenderingIntent::kRelative, &rendering_intent));

    // We didn't have ICC, so all fields should be known.
    if (color_space_ == ColorSpace::kUnknown || tf.IsUnknown()) {
      return JXL_FAILURE(
          "No ICC but cs %u and tf %u%s",
          static_cast<unsigned int>(color_space_),
          tf.IsGamma() ? 0
                       : static_cast<unsigned int>(tf.GetTransferFunction()),
          tf.IsGamma() ? "(gamma)" : "");
    }

    JXL_RETURN_IF_ERROR(CreateICC());
  }

  if (WantICC() && visitor->IsReading()) {
    // Haven't called SetICC() yet, do nothing.
  } else {
    if (ICC().empty()) return JXL_FAILURE("Empty ICC");
  }

  return true;
}


Status ColorEncoding::SetFieldsFromICC() {
  // In case parsing fails, mark the ColorEncoding as invalid.
  SetColorSpace(ColorSpace::kUnknown);
  tf.SetTransferFunction(TransferFunction::kUnknown);

  if (icc_.empty()) return JXL_FAILURE("Empty ICC profile");

#if JPEGXL_ENABLE_SKCMS
  if (icc_.size() < 128) {
    return JXL_FAILURE("ICC file too small");
  }

  skcms_ICCProfile profile;
  JXL_RETURN_IF_ERROR(skcms_Parse(icc_.data(), icc_.size(), &profile));

  // skcms does not return the rendering intent, so get it from the file. It
  // is encoded as big-endian 32-bit integer in bytes 60..63.
  uint32_t rendering_intent32 = icc_[67];
  if (rendering_intent32 > 3 || icc_[64] != 0 || icc_[65] != 0 ||
      icc_[66] != 0) {
    return JXL_FAILURE("Invalid rendering intent %u\n", rendering_intent32);
  }
  // ICC and RenderingIntent have the same values (0..3).
  rendering_intent = static_cast<RenderingIntent>(rendering_intent32);

  if (profile.has_CICP && ApplyCICP(profile.CICP.color_primaries,
                                    profile.CICP.transfer_characteristics,
                                    profile.CICP.matrix_coefficients,
                                    profile.CICP.video_full_range_flag, this)) {
    return true;
  }

  SetColorSpace(ColorSpaceFromProfile(profile));
  cmyk_ = (profile.data_color_space == skcms_Signature_CMYK);

  CIExy wp_unadapted;
  JXL_RETURN_IF_ERROR(UnadaptedWhitePoint(profile, &wp_unadapted));
  JXL_RETURN_IF_ERROR(SetWhitePoint(wp_unadapted));

  // Relies on color_space.
  JXL_RETURN_IF_ERROR(IdentifyPrimaries(profile, wp_unadapted, this));

  // Relies on color_space/white point/primaries being set already.
  DetectTransferFunction(profile, this);
#else  // JPEGXL_ENABLE_SKCMS

  const cmsContext context = GetContext();

  Profile profile;
  JXL_RETURN_IF_ERROR(DecodeProfile(context, icc_, &profile));

  static constexpr size_t kCICPSize = 12;
  static constexpr auto kCICPSignature =
      static_cast<cmsTagSignature>(0x63696370);
  uint8_t cicp_buffer[kCICPSize];
  if (cmsReadRawTag(profile.get(), kCICPSignature, cicp_buffer, kCICPSize) ==
          kCICPSize &&
      ApplyCICP(cicp_buffer[8], cicp_buffer[9], cicp_buffer[10],
                cicp_buffer[11], this)) {
    return true;
  }

  const cmsUInt32Number rendering_intent32 =
      cmsGetHeaderRenderingIntent(profile.get());
  if (rendering_intent32 > 3) {
    return JXL_FAILURE("Invalid rendering intent %u\n", rendering_intent32);
  }
  // ICC and RenderingIntent have the same values (0..3).
  rendering_intent = static_cast<RenderingIntent>(rendering_intent32);

  SetColorSpace(ColorSpaceFromProfile(profile));
  if (cmsGetColorSpace(profile.get()) == cmsSigCmykData) {
    cmyk_ = true;
    return true;
  }

  const cmsCIEXYZ wp_unadapted = UnadaptedWhitePoint(context, profile, *this);
  JXL_RETURN_IF_ERROR(SetWhitePoint(CIExyFromXYZ(wp_unadapted)));

  // Relies on color_space.
  JXL_RETURN_IF_ERROR(IdentifyPrimaries(context, profile, wp_unadapted, this));

  // Relies on color_space/white point/primaries being set already.
  DetectTransferFunction(context, profile, this);

#endif  // JPEGXL_ENABLE_SKCMS

  return true;
}

void ColorEncoding::DecideIfWantICC() {
  PaddedBytes icc_new;
#if JPEGXL_ENABLE_SKCMS
  skcms_ICCProfile profile;
  if (!DecodeProfile(ICC().data(), ICC().size(), &profile)) return;
  if (!MaybeCreateProfile(*this, &icc_new)) return;
#else   // JPEGXL_ENABLE_SKCMS
  const cmsContext context = GetContext();
  Profile profile;
  if (!DecodeProfile(context, ICC(), &profile)) return;
  if (cmsGetColorSpace(profile.get()) == cmsSigCmykData) return;
  if (!MaybeCreateProfile(*this, &icc_new)) return;
#endif  // JPEGXL_ENABLE_SKCMS

  want_icc_ = false;
}


void ConvertInternalToExternalColorEncoding(const ColorEncoding& internal,
                                            JxlColorEncoding* external) {
  external->color_space = static_cast<JxlColorSpace>(internal.GetColorSpace());

  external->white_point = static_cast<JxlWhitePoint>(internal.white_point);

  jxl::CIExy whitepoint = internal.GetWhitePoint();
  external->white_point_xy[0] = whitepoint.x;
  external->white_point_xy[1] = whitepoint.y;

  if (external->color_space == JXL_COLOR_SPACE_RGB ||
      external->color_space == JXL_COLOR_SPACE_UNKNOWN) {
    external->primaries = static_cast<JxlPrimaries>(internal.primaries);
    jxl::PrimariesCIExy primaries = internal.GetPrimaries();
    external->primaries_red_xy[0] = primaries.r.x;
    external->primaries_red_xy[1] = primaries.r.y;
    external->primaries_green_xy[0] = primaries.g.x;
    external->primaries_green_xy[1] = primaries.g.y;
    external->primaries_blue_xy[0] = primaries.b.x;
    external->primaries_blue_xy[1] = primaries.b.y;
  }

  if (internal.tf.IsGamma()) {
    external->transfer_function = JXL_TRANSFER_FUNCTION_GAMMA;
    external->gamma = internal.tf.GetGamma();
  } else {
    external->transfer_function =
        static_cast<JxlTransferFunction>(internal.tf.GetTransferFunction());
    external->gamma = 0;
  }

  external->rendering_intent =
      static_cast<JxlRenderingIntent>(internal.rendering_intent);
}

Status ConvertExternalToInternalColorEncoding(const JxlColorEncoding& external,
                                              ColorEncoding* internal) {
  internal->SetColorSpace(static_cast<ColorSpace>(external.color_space));

  JXL_RETURN_IF_ERROR(ConvertExternalToInternalWhitePoint(
      external.white_point, &internal->white_point));
  if (external.white_point == JXL_WHITE_POINT_CUSTOM) {
    CIExy wp;
    wp.x = external.white_point_xy[0];
    wp.y = external.white_point_xy[1];
    JXL_RETURN_IF_ERROR(internal->SetWhitePoint(wp));
  }

  if (external.color_space == JXL_COLOR_SPACE_RGB ||
      external.color_space == JXL_COLOR_SPACE_UNKNOWN) {
    JXL_RETURN_IF_ERROR(ConvertExternalToInternalPrimaries(
        external.primaries, &internal->primaries));
    if (external.primaries == JXL_PRIMARIES_CUSTOM) {
      PrimariesCIExy primaries;
      primaries.r.x = external.primaries_red_xy[0];
      primaries.r.y = external.primaries_red_xy[1];
      primaries.g.x = external.primaries_green_xy[0];
      primaries.g.y = external.primaries_green_xy[1];
      primaries.b.x = external.primaries_blue_xy[0];
      primaries.b.y = external.primaries_blue_xy[1];
      JXL_RETURN_IF_ERROR(internal->SetPrimaries(primaries));
    }
  }
  CustomTransferFunction tf;
  tf.nonserialized_color_space = internal->GetColorSpace();
  if (external.transfer_function == JXL_TRANSFER_FUNCTION_GAMMA) {
    JXL_RETURN_IF_ERROR(tf.SetGamma(external.gamma));
  } else {
    TransferFunction tf_enum;
    // JXL_TRANSFER_FUNCTION_GAMMA is not handled by this function since there's
    // no internal enum value for it.
    JXL_RETURN_IF_ERROR(ConvertExternalToInternalTransferFunction(
        external.transfer_function, &tf_enum));
    tf.SetTransferFunction(tf_enum);
  }
  internal->tf = tf;

  JXL_RETURN_IF_ERROR(ConvertExternalToInternalRenderingIntent(
      external.rendering_intent, &internal->rendering_intent));

  // The ColorEncoding caches an ICC profile it created earlier that may no
  // longer match the profile with the changed fields, so re-create it.
  if (!(internal->CreateICC())) {
    // This is not an error: for example, it doesn't have ICC profile creation
    // implemented for XYB. This should not be returned as error, since
    // ConvertExternalToInternalColorEncoding still worked correctly, and what
    // matters is that internal->ICC() will not return the wrong profile.
  }

  return true;
}

/* Chromatic adaptation matrices*/
static const float kBradford[9] = {
    0.8951f, 0.2664f, -0.1614f, -0.7502f, 1.7135f,
    0.0367f, 0.0389f, -0.0685f, 1.0296f,
};

static const float kBradfordInv[9] = {
    0.9869929f, -0.1470543f, 0.1599627f, 0.4323053f, 0.5183603f,
    0.0492912f, -0.0085287f, 0.0400428f, 0.9684867f,
};

// Adapts whitepoint x, y to D50
Status AdaptToXYZD50(float wx, float wy, float matrix[9]) {
  if (wx < 0 || wx > 1 || wy <= 0 || wy > 1) {
    // Out of range values can cause division through zero
    // further down with the bradford adaptation too.
    return JXL_FAILURE("Invalid white point");
  }
  float w[3] = {wx / wy, 1.0f, (1.0f - wx - wy) / wy};
  // 1 / tiny float can still overflow
  JXL_RETURN_IF_ERROR(std::isfinite(w[0]) && std::isfinite(w[2]));
  float w50[3] = {0.96422f, 1.0f, 0.82521f};

  float lms[3];
  float lms50[3];

  Mul3x3Vector(kBradford, w, lms);
  Mul3x3Vector(kBradford, w50, lms50);

  if (lms[0] == 0 || lms[1] == 0 || lms[2] == 0) {
    return JXL_FAILURE("Invalid white point");
  }
  float a[9] = {
      //       /----> 0, 1, 2, 3,          /----> 4, 5, 6, 7,          /----> 8,
      lms50[0] / lms[0], 0, 0, 0, lms50[1] / lms[1], 0, 0, 0, lms50[2] / lms[2],
  };
  if (!std::isfinite(a[0]) || !std::isfinite(a[4]) || !std::isfinite(a[8])) {
    return JXL_FAILURE("Invalid white point");
  }

  float b[9];
  Mul3x3Matrix(a, kBradford, b);
  Mul3x3Matrix(kBradfordInv, b, matrix);

  return true;
}

Status PrimariesToXYZ(float rx, float ry, float gx, float gy, float bx,
                      float by, float wx, float wy, float matrix[9]) {
  if (wx < 0 || wx > 1 || wy <= 0 || wy > 1) {
    return JXL_FAILURE("Invalid white point");
  }
  // TODO(lode): also require rx, ry, gx, gy, bx, to be in range 0-1? ICC
  // profiles in theory forbid negative XYZ values, but in practice the ACES P0
  // color space uses a negative y for the blue primary.
  float primaries[9] = {
      rx, gx, bx, ry, gy, by, 1.0f - rx - ry, 1.0f - gx - gy, 1.0f - bx - by};
  float primaries_inv[9];
  memcpy(primaries_inv, primaries, sizeof(float) * 9);
  JXL_RETURN_IF_ERROR(Inv3x3Matrix(primaries_inv));

  float w[3] = {wx / wy, 1.0f, (1.0f - wx - wy) / wy};
  // 1 / tiny float can still overflow
  JXL_RETURN_IF_ERROR(std::isfinite(w[0]) && std::isfinite(w[2]));
  float xyz[3];
  Mul3x3Vector(primaries_inv, w, xyz);

  float a[9] = {
      xyz[0], 0, 0, 0, xyz[1], 0, 0, 0, xyz[2],
  };

  Mul3x3Matrix(primaries, a, matrix);
  return true;
}

Status PrimariesToXYZD50(float rx, float ry, float gx, float gy, float bx,
                         float by, float wx, float wy, float matrix[9]) {
  float toXYZ[9];
  JXL_RETURN_IF_ERROR(PrimariesToXYZ(rx, ry, gx, gy, bx, by, wx, wy, toXYZ));
  float d50[9];
  JXL_RETURN_IF_ERROR(AdaptToXYZD50(wx, wy, d50));

  Mul3x3Matrix(d50, toXYZ, matrix);
  return true;
}

}  // namespace jxl
