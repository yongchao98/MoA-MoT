# Import the kociemba library
import kociemba

# Define the scrambled cube state
scrambled_state = (
    'RBBGYYGYY'
    'BRRWGGOWRR'
    'BRRWGGOBBY'
    'RRYBGGOOOB'
    'OWWYWYWYW'
)

# Solve the cube using the kociemba solver
solution = kociemba.solve(scrambled_state)

# Print the solution in Singmaster notation
print(solution)