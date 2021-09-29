add_subdirectory(new_lib)

# TRIBITS_ADD_LIBRARY(
# core
# HEADERS ${HEADERS}
# SOURCES ${SOURCES}
# )

# TRIBITS_ADD_EXECUTABLE(
# 	Ace_DUMMY
# 	SOURCES main2.cpp #new_lib/src/new_lib.cpp
# 	IMPORTEDLIBS new_lib
# )


TRIBITS_ADD_EXECUTABLE(
	Ace
	SOURCES main.cpp
)

# TRIBITS_ADD_TEST(
# 	Ace
# 	NAME AceTest
# 	ARGS  
# 	NUM_MPI_PROCS 4
# )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(
# 	data_ace_problem
# 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
# 	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersProblemSmall.xml parametersSolver.xml parametersTeko.xml
# 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
# 	DEST_FILES parametersPrec.xml parametersProblem.xml parametersProblemSmall.xml parametersSolver.xml parametersTeko.xml
# 	EXEDEPS Ace
# )


# TRIBITS_COPY_FILES_TO_BINARY_DIR(
# 	data_ace_mesh
#    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
#    SOURCE_FILES tpm_square_h01.mesh
#    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#    DEST_FILES tpm_square_h01.mesh
#    EXEDEPS Ace
# )
