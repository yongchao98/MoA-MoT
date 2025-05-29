# Import the kociemba library
import kociemba

# Define the scrambled state of the cube
scrambled_state = (
    "BBWWYYGOO"
    "GGWWRWWBO"
    "YYOOGBOOR"
    "YRYOGBOWW"
    "GRBWYWRBR"
    "GRBWYWRBR"
)

# Find the solution using the kociemba library
solution = kociemba.solve(scrambled_state)

# Split the solution into individual moves
solution_moves = solution.replace("2", "").replace("'", "").replace(" ", " ")

# Print the solution
print(solution_moves)