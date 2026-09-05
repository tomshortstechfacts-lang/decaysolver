add_executable(decaysolver_cli decaysolver_cli.cpp)
set_target_properties(decaysolver_cli PROPERTIES OUTPUT_NAME decaysolver)
target_link_libraries(decaysolver_cli PRIVATE DecaySolver::decaysolver decaysolver_warnings decaysolver_fp)
# Bibliothèque par défaut : celle du dépôt. Surchargeable par --library.
target_compile_definitions(decaysolver_cli
    PRIVATE DECAYSOLVER_DEFAULT_LIBRARY="${PROJECT_SOURCE_DIR}/data/nuclides_icrp107.csv")
