
TRIBITS_ADD_EXECUTABLE(
	nonconforming
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	nonconforming
    NAME nonconforming
    ARGS
    NUM_MPI_PROCS 4
)

# TRIBITS_COPY_FILES_TO_BINARY_DIR(data_nonconforming
# 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
# 	# SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
# 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
# 	# DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
# 	EXEDEPS nonconforming
# )



# FIXME
find_package(ParMoonolith 0 REQUIRED)
target_link_libraries(problems_nonconforming ParMoonolith::par_moonolith)