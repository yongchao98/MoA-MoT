from kociemba import solve

# Define the scrambled cube state in the kociemba format
# The format is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Each face is represented by 9 characters, starting with the Up face, then Right, Front, Down, Left, Back
# The given cube state is:
# U: OYOBYOBWW
# R: RBBYRYRYO
# F: OGGWGGROW
# D: GYYWWWBRR
# L: GGOYGGROW
# B: YRYYRBBWB

scrambled_state = "OYOBYOBWWRBBYRYRYOGGWGGROWGYYWWWBRRGGOYGGROWYRYYRBBWB"

# Solve the cube
solution = solve(scrambled_state)

# Print the solution
print(solution)