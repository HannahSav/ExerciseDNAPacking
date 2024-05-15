cmake_minimum_required(VERSION 3.10)

# Устанавливаем имя проекта
project(HiCProject)

# Устанавливаем стандарт C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Добавляем исполняемый файл
add_executable(HiCProject main.cpp hicio.cpp counter.cpp)
