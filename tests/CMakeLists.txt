# Catch2 v3 : d'abord une installation système, sinon FetchContent avec un tag figé
# (jamais `main` : un tag flottant rend le build non reproductible).
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    include(FetchContent)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.8.1
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    # Les en-têtes de Catch2 sont traités comme système : nos avertissements stricts (-Werror en
    # CI) portent sur notre code, pas sur celui de la dépendance.
    get_target_property(_catch2_includes Catch2 INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(Catch2 PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_catch2_includes}")
    unset(_catch2_includes)
endif()
include(Catch)

add_subdirectory(unit)
add_subdirectory(regression)
