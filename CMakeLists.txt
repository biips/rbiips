
# configure build
math(EXPR BITS "${CMAKE_SIZEOF_VOID_P} * 8")

if (COMPILER_NAME STREQUAL clang)
    set(RBIIPS_COMPILE_FLAGS "-stdlib=libstdc++ -fPIC -Wno-redeclared-class-member -Wno-c++11-extensions")
    set(RBIIPS_LINK_FLAGS "-stdlib=libstdc++")
  elseif (COMPILER_NAME STREQUAL gcc)
      # TODO : check if g++ version >= 4.7 or 4.8 ?
    set(RBIIPS_COMPILE_FLAGS "-std=c++11 -Wno-unused-local-typedefs")
    set(RBIIPS_LINK_FLAGS ${BIIPS_LINK_FLAGS})

  else()
    message(WARNING "FIXME: Could not determine compiler name")
    set(RBIIPS_COMPILE_FLAGS "-fPIC")
  endif()

elseif (WIN32) # Windows
  set(RBIIPS_COMPILE_FLAGS "-m${BITS}")
  set(RBIIPS_LINK_FLAGS "-m${BITS}")
endif()

# Name of the R package
set (PKG_NAME rbiips)

# add debug flag
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set (R_INSTALL_FLAGS -d)
endif()

if (WIN32)

    if ((${R_ARCH} STREQUAL "x64") AND (BITS EQUAL 32))
        message(FATAL_ERROR "Can not build rbiips for 32bit. R architecture is 64bit. You can use cmake option -DBUILD_RBIIPS=OFF to deactivate rbiips build")
    elseif ((${R_ARCH} STREQUAL "i386") AND (BITS EQUAL 64))
        message(FATAL_ERROR "Can not build rbiips for 64bit. R architecture is 32bit. You can use cmake option -DBUILD_RBIIPS=OFF to deactivate rbiips build")
    endif ()

    # add R CMD INSTALL --build target
    add_custom_target (${PKG_NAME}_build_bin
        DEPENDS ${PKG_NAME}_pkg_dir
        COMMAND ${R_EXECUTABLE} ${R_FLAGS} CMD INSTALL --preclean --build --no-multiarch ${R_INSTALL_FLAGS} ${PKG_BIN_DIR}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building ${PKG_NAME} binary package"
        VERBATIM
    )

elseif (UNIX)

    # add R CMD INSTALL --build target
    add_custom_target (${PKG_NAME}_build_bin
        DEPENDS ${PKG_NAME}_pkg_dir
        COMMAND ${R_EXECUTABLE} ${R_FLAGS} CMD INSTALL --build ${R_INSTALL_FLAGS} ${PKG_BIN_DIR}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building ${PKG_NAME} binary package"
        VERBATIM
    )

endif()

# add R CMD build target
add_custom_target (${PKG_NAME}_build_src
    DEPENDS ${PKG_NAME}_pkg_dir
    COMMAND ${R_EXECUTABLE} ${R_FLAGS} CMD build ${PKG_BIN_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building ${PKG_NAME} source package"
    VERBATIM
)

# add R CMD Rd2pdf target
add_custom_target (${PKG_NAME}_rd2pdf
    DEPENDS ${PKG_NAME}_extract_src
    COMMAND ${R_EXECUTABLE} ${R_FLAGS} CMD Rd2pdf --no-preview --force ${PKG_NAME}_src/${PKG_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  	COMMENT "Generate ${PKG_NAME} documentation PDF output from the Rd sources"
  	VERBATIM
)
