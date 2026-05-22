find_package(Catch2 REQUIRED)

add_executable(test_orbit test_orbital_elements.cpp)

target_link_libraries(test_orbit PRIVATE
    orbit
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(test_orbit)