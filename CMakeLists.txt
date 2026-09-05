cmake_minimum_required(VERSION 3.21)

project(decaysolver
    VERSION 0.1.0
    DESCRIPTION "Résolution vérifiée des équations de Bateman (chaînes de décroissance radioactive)"
    HOMEPAGE_URL "https://github.com/tomshortstechfacts-lang/decaysolver"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options : activées quand decaysolver est le projet principal, désactivées quand il est
# consommé comme sous-projet (add_subdirectory / FetchContent).
option(DECAYSOLVER_BUILD_TESTS "Compiler les tests (Catch2)" ${PROJECT_IS_TOP_LEVEL})
option(DECAYSOLVER_BUILD_CLI "Compiler l'exécutable en ligne de commande" ${PROJECT_IS_TOP_LEVEL})
option(DECAYSOLVER_WARNINGS_AS_ERRORS "Traiter les avertissements comme des erreurs (activé en CI)" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CompilerWarnings)   # cible INTERFACE decaysolver_warnings
include(FloatingPoint)      # cible INTERFACE decaysolver_fp ; refuse -ffast-math
include(GitVersion)         # génère decaysolver/version.hpp

add_subdirectory(src)

if(DECAYSOLVER_BUILD_CLI)
    add_subdirectory(apps)
endif()

if(DECAYSOLVER_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()
