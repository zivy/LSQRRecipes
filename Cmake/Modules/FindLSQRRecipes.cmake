# 
#Find a LSQRRecipes installation or build tree.
#Module is based on the FindITK.cmake module.
#

#search only if the location is not already known.
if(NOT LSQRRecipes_DIR)
  #get the system search path as a list.
  if(UNIX)
    string(REGEX MATCHALL "[^:]+" LSQRRecipes_DIR_SEARCH1 "$ENV{PATH}")
  else()
    string(REGEX REPLACE "\\\\" "/" LSQRRecipes_DIR_SEARCH1 "$ENV{PATH}")
  endif()
  string(REGEX REPLACE "/;" ";" LSQRRecipes_DIR_SEARCH2 ${LSQRRecipes_DIR_SEARCH1})

  #construct a set of paths relative to the system search path.
  set(LSQRRecipes_DIR_SEARCH "")
  foreach(dir ${LSQRRecipes_DIR_SEARCH2})
    set(LSQRRecipes_DIR_SEARCH ${LSQRRecipes_DIR_SEARCH} "${dir}/../lib/LSQRRecipes")
  endforeach(dir)

  #
  #search for an installation or build tree.
  #
  find_path(LSQRRecipes_DIR LSQRRecipesConfig.cmake
            #look for an environment variable LSQRRecipes_DIR
    $ENV{LSQRRecipes_DIR}

             #look in places relative to the system executable search path.
    ${LSQRRecipes_DIR_SEARCH}

             #look in standard UNIX install locations.
    /usr/local/lib/LSQRRecipes
    /usr/lib/LSQRRecipes

           #read from the CMakeSetup registry entries.  It is likely that
           #LSQRRecipes will have been recently built.
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]

            #help the user find it if we cannot.
    DOC "The directory containing LSQRRecipesConfig.cmake. This is either the root of the build tree, or the installation."
  )
endif()

   #if LSQRRecipes was found, load the configuration file to get the rest of the
   #settings.
if(LSQRRecipes_DIR)
  set(LSQRRecipes_FOUND 1)
  include(${LSQRRecipes_DIR}/LSQRRecipesConfig.cmake)

          #set USE_LSQRRecipes_FILE for backward-compatability.
  set(USE_LSQRRecipes_FILE ${LSQRRecipes_USE_FILE})
else(LSQRRecipes_DIR)
  set(LSQRRecipes_FOUND 0)
  if(LSQRRecipes_FIND_REQUIRED)
    message(FATAL_ERROR "Please set LSQRRecipes_DIR to the directory containing LSQRRecipesConfig.cmake. This is either the root of the build tree, or the installation.")
  endif()
endif()
