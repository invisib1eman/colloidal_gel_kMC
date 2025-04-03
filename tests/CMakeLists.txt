# tests/CMakeLists.txt
find_package(GTest REQUIRED)

# Add test executable
add_executable(unit_tests
    test_main.cpp
    test_system.cpp
    # Add other test files
)

target_link_libraries(unit_tests
    PRIVATE
        GTest::GTest
        GTest::Main
        ${PROJECT_NAME}
)

# Register tests
include(GoogleTest)
gtest_discover_tests(unit_tests)