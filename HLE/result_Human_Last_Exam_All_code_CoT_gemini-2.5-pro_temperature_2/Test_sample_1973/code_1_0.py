# First, ensure you have the PuLP library installed.
# You can install it by running the following command in your shell:
# pip install pulp

import pulp

# This plan determines the minimum number of unicorns needed to attack all
# black squares on an 8x8x8 3D chessboard using Integer Linear Programming.

# Step 1: Define the chessboard and identify the black squares.
# An 8x8x8 board has coordinates ranging from 1 to 8.
BOARD_SIZE = 8
all_squares = [(x, y, z) for x in range(1, BOARD_SIZE + 1)
                        for y in range(1, BOARD_SIZE + 1)
                        for z in range(1, BOARD_SIZE + 1)]

# A square (x, y, z) is considered black if the sum of its coordinates is even.
black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 == 0]

# Step 2: Set up the Integer Linear Programming (ILP) problem.
# We aim to find a minimum, so we use LpMinimize.
problem = pulp.LpProblem("Unicorn_Domination", pulp.LpMinimize)

# Step 3: Define the decision variables for the model.
# A binary variable is created for each square.
# place_unicorn[square] = 1 indicates a unicorn is on that square, 0 otherwise.
place_unicorn = pulp.LpVariable.dicts("place_unicorn", all_squares, cat='Binary')

# Step 4: Define the objective function.
# The goal is to minimize the total number of placed unicorns, which is the sum
# of all 'place_unicorn' decision variables.
problem += pulp.lpSum([place_unicorn[s] for s in all_squares])

# Step 5: Define the constraints for the model.
# Every black square must be attacked by at least one unicorn.
# A unicorn at square 'i' attacks square 'j' if they are on a space diagonal,
# which is true if abs(i.x - j.x) == abs(i.y - j.y) == abs(i.z - j.z).

print("Setting up the ILP constraints. This might take a few moments...")
for j in black_squares:
    # For each black square 'j', we identify all squares 'i' from which a unicorn can attack it.
    attacking_squares_vars = [
        place_unicorn[i] for i in all_squares
        if abs(i[0] - j[0]) == abs(i[1] - j[1]) == abs(i[2] - j[2])
    ]
    
    # We add the constraint that the sum of unicorns on these attacking squares must be >= 1.
    problem += pulp.lpSum(attacking_squares_vars) >= 1

# Step 6: Solve the ILP problem.
print("Solving the ILP problem. This may take some time depending on your system...")
# We use the CBC solver that is included with PuLP.
# The 'msg=False' argument suppresses the solver's log output.
problem.solve(pulp.PULP_CBC_CMD(msg=False))

# Step 7: Print the result of the optimization.
if problem.status == pulp.LpStatusOptimal:
    # Extract the optimal value for the objective function.
    min_unicorns = int(pulp.value(problem.objective))
    
    print("\n-------------------")
    print("      RESULT       ")
    print("-------------------")
    print(f"The ILP problem was solved to optimality.")
    print(f"The minimum number of unicorns required is: {min_unicorns}")

    # As requested, output the final answer as an equation where each '1' represents one unicorn.
    print("\nThe final equation is:")
    sum_string = " + ".join(["1"] * min_unicorns)
    print(f"{sum_string} = {min_unicorns}")
else:
    print(f"\nThe optimal solution could not be found. Status: {pulp.LpStatus[problem.status]}")
