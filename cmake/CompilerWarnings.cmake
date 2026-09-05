# Cible INTERFACE portant les avertissements stricts du projet.
# À lier en PRIVATE sur chaque cible compilée du projet (jamais exportée aux consommateurs).

add_library(decaysolver_warnings INTERFACE)

if(MSVC)
    target_compile_options(decaysolver_warnings INTERFACE
        /W4
        /permissive-
        /utf-8
        /w14242 # conversion avec perte possible de données
        /w14254 # conversion d'un champ de bits avec perte
        /w14263 # fonction membre qui ne redéfinit aucune fonction virtuelle
        /w14265 # classe avec fonctions virtuelles sans destructeur virtuel
        /w14287 # comparaison non signée / négatif
        /we4289 # variable de boucle utilisée hors de la boucle
        /w14296 # expression toujours vraie/fausse
        /w14311 # troncature de pointeur
        /w14545 /w14546 /w14547 /w14549 /w14555 # expressions sans effet
        /w14619 # pragma warning sur un numéro inexistant
        /w14640 # objet statique local non thread-safe
        /w14826 # extension de signe d'une conversion
        /w14905 /w14906 # cast de chaîne vers LPSTR/LPWSTR
        /w14928 # conversion de copie illégale
    )
    if(DECAYSOLVER_WARNINGS_AS_ERRORS)
        target_compile_options(decaysolver_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(decaysolver_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        -Wsign-conversion
        -Wfloat-equal        # interdit `==` entre flottants : les comparaisons passent par des tolérances
        -Wshadow
        -Wold-style-cast
        -Wcast-align
        -Wnon-virtual-dtor
        -Wdouble-promotion
        -Wformat=2
        -Wimplicit-fallthrough
        -Wnull-dereference
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(decaysolver_warnings INTERFACE
            -Wduplicated-cond
            -Wduplicated-branches
            -Wlogical-op
            -Wuseless-cast
        )
    endif()
    if(DECAYSOLVER_WARNINGS_AS_ERRORS)
        target_compile_options(decaysolver_warnings INTERFACE -Werror)
    endif()
endif()
