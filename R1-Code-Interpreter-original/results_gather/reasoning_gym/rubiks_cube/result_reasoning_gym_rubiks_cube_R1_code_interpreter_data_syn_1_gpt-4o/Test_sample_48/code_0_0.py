# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'YYOYYWGG',  # Up face
    'RRWOOYYG',  # Left face
    'WRRGGWOO',  # Front face
    'YRGBBWOB',  # Right face
    'WBRWGYBR',  # Back face
    'BWOOWOYR'   # Down face
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' ').replace('\'', ' \'').replace('  ', ' ')

# Print the solution
print(formatted_solution)