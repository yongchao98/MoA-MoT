import kociemba

# Define the cube state
cube_state = (
    "OOB"
    "GYY"
    "GGG"
    "YRRWWWOO"
    "YRWO"
    "GGYO"
    "WRBB"
    "GRWO"
    "GGYB"
    "BYY"
    "BOO"
    "BWW"
    "RRR"
)

# Solve the cube
solution = kociemba.solve(cube_state)

# Convert the solution to the required format
solution_steps = solution.replace("2", " ").replace("'", " ").split()
formatted_solution = ' '.join(solution_steps)

print(formatted_solution)