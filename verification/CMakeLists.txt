# Exécutables de vérification : ils produisent les tableaux du rapport et rendent un code de
# retour exploitable par la CI. Ils partagent l'assistant de chaîne des tests unitaires.

add_executable(decaysolver_convergence V2_order_of_accuracy/convergence.cpp)
target_include_directories(decaysolver_convergence PRIVATE "${PROJECT_SOURCE_DIR}/tests/unit")
target_link_libraries(decaysolver_convergence
    PRIVATE DecaySolver::decaysolver decaysolver_warnings decaysolver_fp)
target_compile_definitions(decaysolver_convergence PRIVATE
    DECAYSOLVER_CONVERGENCE_ORACLE="${PROJECT_SOURCE_DIR}/verification/V2_order_of_accuracy/oracle_lambda123.csv")

# Le critère d'acceptation des ordres de convergence est aussi un test ctest : une régression de
# la numérique qui casserait un ordre fait échouer la suite de tests.
add_test(NAME verification.V2.convergence_orders COMMAND decaysolver_convergence)
set_tests_properties(verification.V2.convergence_orders PROPERTIES LABELS "verification")
