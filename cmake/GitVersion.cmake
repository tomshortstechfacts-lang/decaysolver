# Génère include/decaysolver/version.hpp (dans l'arborescence de build) à partir de la version
# du projet et de l'état du dépôt git, pour l'en-tête de provenance des sorties.
#
# Limite connue : le SHA est capturé à la CONFIGURATION, pas à chaque build. Après un commit,
# relancer la configuration (cmake --preset ...) pour que les sorties portent le bon SHA.
# L'en-tête de provenance le rappelle explicitement.

find_package(Git QUIET)

set(DECAYSOLVER_GIT_SHA "inconnu")
set(DECAYSOLVER_GIT_DIRTY "false")

if(Git_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE _sha
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _rc
        ERROR_QUIET)
    if(_rc EQUAL 0)
        set(DECAYSOLVER_GIT_SHA "${_sha}")
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" diff --quiet HEAD --
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            RESULT_VARIABLE _dirty_rc
            ERROR_QUIET)
        if(NOT _dirty_rc EQUAL 0)
            set(DECAYSOLVER_GIT_DIRTY "true")
        endif()
    endif()
    unset(_sha)
    unset(_rc)
    unset(_dirty_rc)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(DECAYSOLVER_BUILD_TYPE "non spécifié (générateur multi-configuration)")
else()
    set(DECAYSOLVER_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
endif()

configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/version.hpp.in"
    "${PROJECT_BINARY_DIR}/generated/decaysolver/version.hpp"
    @ONLY)
