cmake_minimum_required(VERSION 3.5)

project(asu-p3)

add_executable(asu-p3
    src/Fluid2.h
    src/Fluid2Exercise.cpp
    src/FluidVisualizer2.cpp
    src/FluidVisualizer2.h
    src/Grid2.h
    src/main.cpp
    src/Scene.cpp
    src/Scene.h
    src/Containers/Array2.h
    src/Containers/Index2.cpp
    src/Containers/Index2.h
    src/Math/AABox2.h
    src/Math/Vector2.cpp
    src/Math/Vector2.h
    src/Math/Vector3.cpp
    src/Math/Vector3.h
    src/Numeric/SparseMatrix.h
    src/Numeric/PCGSolver.h
    src/Std/CMath.h
    src/Std/CStdInt.h
)

find_package(OpenGL REQUIRED)
find_package(FreeGLUT CONFIG REQUIRED)

target_include_directories(asu-p3
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/src
)

target_link_libraries(asu-p3
    ${OPENGL_LIBRARIES}
    $<IF:$<TARGET_EXISTS:FreeGLUT::freeglut>,FreeGLUT::freeglut,FreeGLUT::freeglut_static>
)
