# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'GGGYYYOOO',  # Top face
    'OBBRRYGGR',  # Middle layer
    'WWWRRYBBB',  # Bottom face
    'RRRWWWBBB',  # Front face
    'GGGWWWBBB',  # Back face
    'OOOYYYRRR'   # Left face
)

# Convert the scrambled state to a string format required by the solver
scrambled_string = ''.join(scrambled_state)

# Solve the cube
solution = utils.solve(scrambled_string, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))