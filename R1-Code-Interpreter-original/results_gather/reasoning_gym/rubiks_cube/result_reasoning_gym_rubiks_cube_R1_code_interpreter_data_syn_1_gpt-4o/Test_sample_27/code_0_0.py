# Import the kociemba library
import kociemba

# Define the scrambled state of the cube
scrambled_state = (
    'RRBYYYYYY'
    'OOOBBBRRW'
    'RRWRGGYOO'
    'BBWGGYGGG'
    'GOOWWWWWW'
)

# Use the kociemba library to find the solution
solution = kociemba.solve(scrambled_state)

# Print the solution
print(solution)