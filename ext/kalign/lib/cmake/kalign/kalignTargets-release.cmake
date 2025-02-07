#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "kalign::kalign" for configuration "Release"
set_property(TARGET kalign::kalign APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(kalign::kalign PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libkalign.a"
  )

list(APPEND _cmake_import_check_targets kalign::kalign )
list(APPEND _cmake_import_check_files_for_kalign::kalign "${_IMPORT_PREFIX}/lib/libkalign.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
