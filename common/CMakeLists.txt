project(common)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Создаем библиотеку
add_library(common STATIC
    src/matrix.cpp
)

# Устанавливаем пути для include файлов
target_include_directories(common PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Устанавливаем include файлы для цели install
install(DIRECTORY include/ DESTINATION include)
install(TARGETS common DESTINATION lib)