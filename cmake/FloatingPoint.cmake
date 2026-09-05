# Sémantique flottante du projet.
#
# Tout le raisonnement sur l'erreur d'arrondi (tolérances des tests, analyse de l'annulation
# catastrophique dans la formule de Bateman) suppose une arithmétique IEEE-754 stricte,
# opération par opération, dans l'ordre écrit dans le source. Deux options de compilation la
# cassent :
#   - -ffast-math (et ses composantes -funsafe-math-optimizations, -ffinite-math-only...) :
#     autorise la réassociation, suppose l'absence de NaN/inf, et rend les résultats
#     dépendants du compilateur et du niveau d'optimisation ;
#   - la contraction en FMA (a*b+c calculé avec un seul arrondi) : change l'arrondi des
#     sommes de produits, donc les résultats, entre machines avec et sans FMA.
# On refuse la première et on désactive la seconde.

add_library(decaysolver_fp INTERFACE)

if(MSVC)
    target_compile_options(decaysolver_fp INTERFACE /fp:precise)
else()
    target_compile_options(decaysolver_fp INTERFACE -ffp-contract=off)
endif()

foreach(_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_DEBUG
             CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_MINSIZEREL)
    if("${${_var}}" MATCHES "-ffast-math|-Ofast|/fp:fast|-funsafe-math-optimizations")
        message(FATAL_ERROR
            "${_var} contient une option qui casse la sémantique IEEE-754 (${${_var}}). "
            "Ce projet la refuse : voir README, section « Hypothèses ».")
    endif()
endforeach()
unset(_var)
