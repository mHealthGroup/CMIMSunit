# add_subdirectory(signal)
# add_subdirectory(stats)
# add_subdirectory(trapz)

# add_executable(
#     mimsunit
#     mims_aggregate.c mims_aggregate.h
#     mims_combine_axes.c mims_combine_axes.h
#     mims_extrapolate.c mims_extrapolate.h
#     mims_filter.c mims_filter.h
#     mims_helper.c mims_helper.h
#     mims_unit.c mims_unit.h
# )

# target_include_directories(mimsunit PUBLIC
#                           "${PROJECT_BINARY_DIR}/src"
#                           "${PROJECT_SOURCE_DIR}/src/signal"
#                           "${PROJECT_SOURCE_DIR}/src/stats"
#                           "${PROJECT_SOURCE_DIR}/src/trapz"
#                           )

# cmake_minimum_required(VERSION 3.5)
# # Project name
# project(CMIMS)

add_executable(
    mimsunit
    mims_aggregate.c
    mims_aggregate.h
    mims_combine_axes.c
    mims_combine_axes.h
    mims_extrapolate.c
    mims_extrapolate.h
    mims_filter.c
    mims_filter.h
    mims_helper.c
    mims_helper.h
    mims_unit.c
    mims_unit.h
    )
    
add_subdirectory(signal)
add_subdirectory(stats)
add_subdirectory(trapz)

target_link_libraries(mimsunit PUBLIC signal stats trapz)

target_include_directories(mimsunit PUBLIC
                          "${PROJECT_SOURCE_DIR}"
                          "${PROJECT_SOURCE_DIR}/signal"
                          "${PROJECT_SOURCE_DIR}/stats"
                          "${PROJECT_SOURCE_DIR}/trapz"
                          )