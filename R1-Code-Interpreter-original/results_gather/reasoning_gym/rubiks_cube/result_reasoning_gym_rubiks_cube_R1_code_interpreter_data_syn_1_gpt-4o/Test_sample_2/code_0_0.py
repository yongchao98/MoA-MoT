from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'Y', 'G', 'G', 'G', 'Y', 'Y', 'G', 'Y', 'Y',  # Top face (U)
    'B', 'R', 'R', 'W', 'R', 'R', 'G', 'G', 'W',  # Left face (L)
    'B', 'R', 'R', 'W', 'G', 'W', 'O', 'O', 'O',  # Front face (F)
    'Y', 'Y', 'Y', 'R', 'G', 'R', 'W', 'W', 'W',  # Right face (R)
    'O', 'B', 'O', 'Y', 'B', 'O', 'O', 'B', 'G',  # Back face (B)
    'B', 'W', 'B', 'B', 'W', 'B', 'O', 'R', 'B'   # Bottom face (D)
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

# Print the solution
print(singmaster_solution)