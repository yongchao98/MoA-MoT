from rubik_solver import utils

# Define the cube configuration
cube_configuration = (
    'YGO'
    'RYR'
    'RBB'
    'OWWGYYO'
    'RRWBGGYO'
    'RRWBYRBO'
    'ROY'
    'YWG'
    'YOG'
)

# Solve the cube
solution = utils.solve(cube_configuration, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))