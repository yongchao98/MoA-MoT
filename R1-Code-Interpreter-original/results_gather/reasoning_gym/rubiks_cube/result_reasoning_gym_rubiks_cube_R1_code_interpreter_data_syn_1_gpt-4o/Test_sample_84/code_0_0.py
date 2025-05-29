# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'YYRYYOBBYBBWRRGOWYBOO'
    'BRWGGROWBBOBRWGGROYYG'
    'WWOOOGGWGRWRR'
)

# Solve the cube using the Kociemba algorithm
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' ').replace("'", "' ")

# Print the solution
print(formatted_solution)