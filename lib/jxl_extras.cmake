# Copyright (c) the JPEG XL Project Authors. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

include(compatibility.cmake)
include(jxl_lists.cmake)

# Object library for those parts of extras that do not depend on jxl internals,
# jpegli, or external codec libraries. We will create two versions of these
# object files, one with and one without -DJPEGXL_ENABLE_* definitions for
# external codec support.
list(APPEND JPEGXL_EXTRAS_CORE_SOURCES
  "${JPEGXL_INTERNAL_EXTRAS_SOURCES}"
  "${JPEGXL_INTERNAL_CODEC_JXL_SOURCES}"
  "${JPEGXL_INTERNAL_CODEC_PGX_SOURCES}"
  "${JPEGXL_INTERNAL_CODEC_PNM_SOURCES}"
  "${JPEGXL_INTERNAL_CODEC_NPY_SOURCES}"
)
foreach(LIB jxl_extras_core-obj jxl_extras_core_nocodec-obj)
  add_library("${LIB}" OBJECT "${JPEGXL_EXTRAS_CORE_SOURCES}")
  list(APPEND JXL_EXTRAS_OBJECT_LIBRARIES "${LIB}")
endforeach()
list(APPEND JXL_EXTRAS_OBJECTS $<TARGET_OBJECTS:jxl_extras_core-obj>)

# Object library for those parts of extras that depend on jxl internals, with
# and without codec support.
foreach(LIB jxl_extras_internal-obj jxl_extras_internal_nocodec-obj)
  add_library("${LIB}" OBJECT
    "${JPEGXL_INTERNAL_EXTRAS_FOR_TOOLS_SOURCES}"
  )
  list(APPEND JXL_EXTRAS_OBJECT_LIBRARIES "${LIB}")
endforeach()
list(APPEND JXL_EXTRAS_OBJECTS $<TARGET_OBJECTS:jxl_extras_internal-obj>)

# Object library for encoders/decoders that depend on external codec libraries.
# The actual sources, dependencies and compile definitions will be populated
# later based on what codecs we find.
add_library(jxl_extras_codec-obj OBJECT)
list(APPEND JXL_EXTRAS_OBJECT_LIBRARIES jxl_extras_codec-obj)
list(APPEND JXL_EXTRAS_OBJECTS $<TARGET_OBJECTS:jxl_extras_codec-obj>)
set(JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES)
set(JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS)

find_package(GIF 5.1)
if(GIF_FOUND)
  target_sources(jxl_extras_codec-obj PRIVATE
    extras/dec/gif.cc
    extras/dec/gif.h
  )
  target_include_directories(jxl_extras_codec-obj PRIVATE "${GIF_INCLUDE_DIRS}")
  list(APPEND JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES ${GIF_LIBRARIES})
  list(APPEND JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS -DJPEGXL_ENABLE_GIF=1)
  if(JPEGXL_DEP_LICENSE_DIR)
    configure_file("${JPEGXL_DEP_LICENSE_DIR}/libgif-dev/copyright"
                   ${PROJECT_BINARY_DIR}/LICENSE.libgif COPYONLY)
  endif()  # JPEGXL_DEP_LICENSE_DIR
endif()

find_package(JPEG)
if(JPEG_FOUND)
  target_sources(jxl_extras_codec-obj PRIVATE
    "${JPEGXL_INTERNAL_CODEC_JPG_SOURCES}"
  )
  target_include_directories(jxl_extras_codec-obj PRIVATE
    "${JPEG_INCLUDE_DIRS}"
  )
  list(APPEND JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES ${JPEG_LIBRARIES})
  list(APPEND JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS -DJPEGXL_ENABLE_JPEG=1)
  if(JPEGXL_DEP_LICENSE_DIR)
    configure_file("${JPEGXL_DEP_LICENSE_DIR}/libjpeg-dev/copyright"
                   ${PROJECT_BINARY_DIR}/LICENSE.libjpeg COPYONLY)
  endif()  # JPEGXL_DEP_LICENSE_DIR
endif()

if(JPEG_FOUND AND JPEGXL_ENABLE_JPEGLI)
  add_library(jxl_extras_jpegli-obj OBJECT
    "${JPEGXL_INTERNAL_CODEC_JPEGLI_SOURCES}"
  )
  target_compile_definitions(jxl_extras_jpegli-obj PRIVATE
    -DJPEGXL_ENABLE_JPEGLI=1
  )
  target_include_directories(jxl_extras_jpegli-obj PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/include/jpegli"
  )
  list(APPEND JXL_EXTRAS_OBJECT_LIBRARIES jxl_extras_jpegli-obj)
  list(APPEND JXL_EXTRAS_OBJECTS $<TARGET_OBJECTS:jxl_extras_jpegli-obj>)
endif()

if(NOT JPEGXL_BUNDLE_LIBPNG)
  find_package(PNG)
endif()
if(PNG_FOUND)
  target_sources(jxl_extras_codec-obj PRIVATE
    "${JPEGXL_INTERNAL_CODEC_APNG_SOURCES}"
  )
  target_include_directories(jxl_extras_codec-obj PRIVATE "${PNG_INCLUDE_DIRS}")
  list(APPEND JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES ${PNG_LIBRARIES})
  list(APPEND JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS -DJPEGXL_ENABLE_APNG=1)
  configure_file(extras/LICENSE.apngdis
                 ${PROJECT_BINARY_DIR}/LICENSE.apngdis COPYONLY)
endif()

if (JPEGXL_ENABLE_OPENEXR)
pkg_check_modules(OpenEXR IMPORTED_TARGET OpenEXR)
if (OpenEXR_FOUND)
  target_sources(jxl_extras_codec-obj PRIVATE
    "${JPEGXL_INTERNAL_CODEC_EXR_SOURCES}"
  )
  list(APPEND JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS -DJPEGXL_ENABLE_EXR=1)
  target_include_directories(jxl_extras_codec-obj PRIVATE
    "${OpenEXR_INCLUDE_DIRS}"
  )
  list(APPEND JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES PkgConfig::OpenEXR)
  if(JPEGXL_DEP_LICENSE_DIR)
    configure_file("${JPEGXL_DEP_LICENSE_DIR}/libopenexr-dev/copyright"
                   ${PROJECT_BINARY_DIR}/LICENSE.libopenexr COPYONLY)
  endif()  # JPEGXL_DEP_LICENSE_DIR
  # OpenEXR generates exceptions, so we need exception support to catch them.
  # Actually those flags counteract the ones set in JPEGXL_INTERNAL_FLAGS.
  if (NOT WIN32)
    set_source_files_properties(
      extras/dec/exr.cc extras/enc/exr.cc PROPERTIES COMPILE_FLAGS -fexceptions)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set_source_files_properties(
	extras/dec/exr.cc extras/enc/exr.cc PROPERTIES COMPILE_FLAGS
	-fcxx-exceptions)
    endif()
  endif()
endif() # OpenEXR_FOUND
endif() # JPEGXL_ENABLE_OPENEXR

# Common settings for the object libraries.
foreach(LIB ${JXL_EXTRAS_OBJECT_LIBRARIES})
  target_compile_options("${LIB}" PRIVATE "${JPEGXL_INTERNAL_FLAGS}")
  target_compile_definitions("${LIB}" PRIVATE -DJXL_EXPORT=)
  set_property(TARGET "${LIB}" PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_include_directories("${LIB}" PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
    ${JXL_HWY_INCLUDE_DIRS}
  )
endforeach()

foreach(LIB jxl_extras_core-obj jxl_extras_internal-obj jxl_extras_codec-obj)
  target_compile_definitions("${LIB}" PUBLIC
    ${JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS}
  )
endforeach()

# Define an extras library that does not have the image codecs, only the core
# extras code. This is needed for some of the fuzzers.
add_library(jxl_extras_nocodec-static STATIC EXCLUDE_FROM_ALL
  $<TARGET_OBJECTS:jxl_extras_core_nocodec-obj>
  $<TARGET_OBJECTS:jxl_extras_internal_nocodec-obj>
)
target_link_libraries(jxl_extras_nocodec-static PUBLIC
  jxl-static
  jxl_threads-static
)

# We only define a static library jxl_extras since it uses internal parts of
# jxl library which are not accessible from outside the library in the
# shared library case.
add_library(jxl_extras-static STATIC EXCLUDE_FROM_ALL ${JXL_EXTRAS_OBJECTS})
target_compile_definitions(jxl_extras-static PUBLIC
  ${JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS}
)
target_link_libraries(jxl_extras-static PUBLIC
  ${JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES}
  jxl-static
  jxl_threads-static
)
if (JPEGXL_ENABLE_SJPEG)
  target_compile_definitions(jxl_extras-static PUBLIC -DJPEGXL_ENABLE_SJPEG=1)
  target_link_libraries(jxl_extras-static PRIVATE sjpeg)
endif ()
if(JPEG_FOUND AND JPEGXL_ENABLE_JPEGLI)
  target_compile_definitions(jxl_extras-static PUBLIC -DJPEGXL_ENABLE_JPEGLI=1)
  target_link_libraries(jxl_extras-static PRIVATE jpegli-static)
endif()

### Static library that does not depend on internal parts of jxl library.
add_library(jxl_extras_codec-static STATIC
  $<TARGET_OBJECTS:jxl_extras_core-obj>
  $<TARGET_OBJECTS:jxl_extras_codec-obj>
)
target_compile_definitions(jxl_extras_codec-static PUBLIC
  ${JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS}
)
target_link_libraries(jxl_extras_codec-static PRIVATE
  ${JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES}
  jxl
)

### Shared library that does not depend on internal parts of jxl library.
### Used by cjxl and djxl binaries.
if (BUILD_SHARED_LIBS)
add_library(jxl_extras_codec SHARED
  $<TARGET_OBJECTS:jxl_extras_codec-obj>
  $<TARGET_OBJECTS:jxl_extras_core-obj>
)
target_compile_definitions(jxl_extras_codec PUBLIC
  ${JXL_EXTRAS_CODEC_PUBLIC_DEFINITIONS}
)
target_link_libraries(jxl_extras_codec PRIVATE
  ${JXL_EXTRAS_CODEC_INTERNAL_LIBRARIES}
  jxl
)
set_target_properties(jxl_extras_codec PROPERTIES
  VERSION ${JPEGXL_LIBRARY_VERSION}
  SOVERSION ${JPEGXL_LIBRARY_SOVERSION}
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
install(TARGETS jxl_extras_codec
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
else()
add_library(jxl_extras_codec ALIAS jxl_extras_codec-static)
endif()  # BUILD_SHARED_LIBS

#TODO(lode): don't depend on CMS for the core library
if (JPEGXL_ENABLE_SKCMS)
  target_include_directories(jxl_extras_core_nocodec-obj PRIVATE
    $<TARGET_PROPERTY:skcms,INCLUDE_DIRECTORIES>
  )
  target_include_directories(jxl_extras_core-obj PRIVATE
    $<TARGET_PROPERTY:skcms,INCLUDE_DIRECTORIES>
  )
else ()
  target_include_directories(jxl_extras_core_nocodec-obj PRIVATE
    $<TARGET_PROPERTY:lcms2,INCLUDE_DIRECTORIES>
  )
endif ()
