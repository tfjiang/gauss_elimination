project(test)
file(GLOB source *.cpp)
file(GLOB head *.h)
add_executable(test ${source} ${header})