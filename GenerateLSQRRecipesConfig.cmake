#
# Generate the LSQRRecipesConfig.cmake file for the build tree and installation.
#

################################################################################
#
# Build tree.
#

set(LSQRRecipes_USE_FILE ${LSQRRecipes_BINARY_DIR}/UseLSQRRecipes.cmake)
set(LSQRRecipes_LIBRARY_DIRS_CONFIG ${LSQRRecipes_BINARY_DIR})
set(LSQRRecipes_INCLUDE_DIRS_CONFIG ${LSQRRecipes_INCLUDE_DIRS})
configure_file(${LSQRRecipes_SOURCE_DIR}/LSQRRecipesConfig.cmake.in
               ${LSQRRecipes_BINARY_DIR}/LSQRRecipesConfig.cmake @ONLY)

################################################################################
#
# Install tree.
#

set(LSQRRecipes_USE_FILE ${LSQRRecipes_INSTALL_LIB_DIR}/UseLSQRRecipes.cmake)
set(LSQRRecipes_LIBRARY_DIRS_CONFIG ${LSQRRecipes_INSTALL_LIB_DIR})
set(LSQRRecipes_INCLUDE_DIRS_CONFIG ${LSQRRecipes_INSTALL_INCLUDE_DIR})

#place the LSQRRecipesConfig.cmake for the installation version in a temporary
#location
configure_file(${LSQRRecipes_SOURCE_DIR}/LSQRRecipesConfig.cmake.in
               ${PRE_INSTALL_DIR}/LSQRRecipesConfig.cmake @ONLY)
