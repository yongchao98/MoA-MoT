# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'GGWGGYWYWWWRGBBOOOGYR'
    'YRRRGGGOOWRBBYYOGRYWR'
    'YBORWOR'
)

# Solve the cube using the Kociemba algorithm
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))