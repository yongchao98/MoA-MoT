# Import the kociemba library
import kociemba

# Define the scrambled state of the cube
scrambled_state = (
    'RRYYGYYGOO'
    'BYYRGGWOGY'
    'WRRWGGWOYB'
    'WRRWGGWOYB'
    'BRR'
    'BWW'
    'OOO'
)

# Use the kociemba solver to find the solution
solution = kociemba.solve(scrambled_state)

# Print the solution in Singmaster notation
print(solution)