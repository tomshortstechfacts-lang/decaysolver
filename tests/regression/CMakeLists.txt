add_executable(decaysolver_regression_tests test_regression.cpp)
target_link_libraries(decaysolver_regression_tests
    PRIVATE DecaySolver::decaysolver Catch2::Catch2WithMain decaysolver_warnings decaysolver_fp)
target_compile_definitions(decaysolver_regression_tests PRIVATE
    DECAYSOLVER_DATA_DIR="${PROJECT_SOURCE_DIR}/data"
    DECAYSOLVER_REGRESSION_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
catch_discover_tests(decaysolver_regression_tests)

# Outil de régénération de la référence, jamais lancé automatiquement.
add_executable(decaysolver_regenerate_reference regenerate_reference.cpp)
target_link_libraries(decaysolver_regenerate_reference
    PRIVATE DecaySolver::decaysolver decaysolver_warnings decaysolver_fp)
target_compile_definitions(decaysolver_regenerate_reference PRIVATE
    DECAYSOLVER_DEFAULT_LIBRARY="${PROJECT_SOURCE_DIR}/data/nuclides_icrp107.csv")
