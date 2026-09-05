add_executable(decaysolver_unit_tests
    test_units.cpp
    test_provenance.cpp
    test_nuclide.cpp
    test_nuclide_library.cpp
    test_decay_system.cpp
    test_bateman.cpp
    test_integrator.cpp
)
target_link_libraries(decaysolver_unit_tests
    PRIVATE DecaySolver::decaysolver Catch2::Catch2WithMain decaysolver_warnings decaysolver_fp)

# Les tests de données lisent data/ depuis l'arborescence source (jamais copié dans le build).
target_compile_definitions(decaysolver_unit_tests
    PRIVATE DECAYSOLVER_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

catch_discover_tests(decaysolver_unit_tests)
