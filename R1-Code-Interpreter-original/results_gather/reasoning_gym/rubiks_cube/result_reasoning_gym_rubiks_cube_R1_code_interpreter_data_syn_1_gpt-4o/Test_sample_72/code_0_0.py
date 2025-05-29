# Importing a Rubik's Cube solver library
from rubik_solver import utils

# The scrambled state of the cube in a format that the solver can understand
# The format is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Each face is represented by 9 characters, starting with the Up face, then Right, Front, Down, Left, Back
# The given state needs to be converted to this format

# Given state:
# U: B G G B Y Y O Y Y
# R: R R Y G R R G G W
# F: O O W R R B Y G W
# D: O O O Y B W R W B
# L: R R B Y G W O O O
# B: Y B W G W B G W B

# Converting the given state to the required format
scrambled_state = "BGGBYYOYY" + "RRYGRRGGW" + "OOWRRBYGW" + "OOOYBWRWB" + "RRBYGWO" + "YBWGWB"

# Solving the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Converting the solution to the required format (no combined steps)
solution_steps = ' '.join(move.replace('2', ' ') for move in solution).replace('  ', ' ')

# Output the solution
print(solution_steps)