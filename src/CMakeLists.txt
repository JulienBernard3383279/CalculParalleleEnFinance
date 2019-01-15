add_library(mc-pricer
    SHARED
    Parser.cpp
    BasketOption.cpp
    AsianOption.cpp
    PerformanceOption.cpp
    BlackScholesModel.cpp
    MonteCarlo.cpp
    )

install(TARGETS mc-pricer
    DESTINATION lib)

file(
    GLOB
    headers
    *.hpp
    )

install(FILES ${headers}
    DESTINATION include/${CMAKE_PROJECT_NAME})

add_executable(pricer
    Pricer.cpp
)

# Extra linking for the project
target_link_libraries(pricer mc-pricer ${LIBS} ${MPI_LIBS})
