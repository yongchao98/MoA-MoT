# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'YRGGYGRYY'
    'BRBWBBROOY'
    'RRYGGBWOOYBB'
    'RGGOOGWWBWO'
    'WWRWWRYBO'
)

# Solve the cube using the library
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

# Print the solution
print(singmaster_solution)