# labx/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(labs)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Добавляем исполняемый файл
add_executable(labs src/main.cpp)

# Связываем исполняемый файл с библиотекой common
target_link_libraries(labs PRIVATE common)
