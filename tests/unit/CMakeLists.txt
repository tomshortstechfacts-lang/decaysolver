add_executable(decaysolver_unit_tests
    test_units.cpp
    test_provenance.cpp
)
target_link_libraries(decaysolver_unit_tests
    PRIVATE DecaySolver::decaysolver Catch2::Catch2WithMain decaysolver_warnings decaysolver_fp)

catch_discover_tests(decaysolver_unit_tests)
