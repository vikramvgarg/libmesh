#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex5

example_dir=examples/adjoints/$example_name

# Save previous timestep solution in memory
run_example "$example_name" solutionhistorytype=memory_solution_history
# Save previous timestep solution on disk
run_example "$example_name" solutionhistorytype=file_solution_history
