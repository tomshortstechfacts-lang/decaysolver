add_library(decaysolver
    units.cpp
    provenance.cpp
)
add_library(DecaySolver::decaysolver ALIAS decaysolver)

target_include_directories(decaysolver
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(decaysolver PUBLIC cxx_std_20)
target_link_libraries(decaysolver PRIVATE decaysolver_warnings decaysolver_fp)
