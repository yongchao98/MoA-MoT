import pulp

# This script might take a minute or two to run as it solves a complex optimization problem.

# Step 1: Explain the plan
print("Modeling the problem as an Integer Linear Program (Set Cover).")
print("Objective: Minimize the number of unicorns.")
print("Constraints: Every black square must be attacked by at least one unicorn.")
print("Solving for an 8x8x8 board...\n")

# Step 2: Define board parameters
N = 8
cells = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

# Step 3: Identify the black squares to be covered
# A square is black if the sum of its coordinates is odd.
black_squares = [c for c in cells if (c[0] + c[1] + c[2]) % 2 != 0]

# Step 4: Create the ILP problem
prob = pulp.LpProblem("Unicorn_Covering", pulp.LpMinimize)

# Step 5: Create the decision variables
# u_c = 1 if a unicorn is at cell c, 0 otherwise
unicorns = pulp.LpVariable.dicts("unicorn", cells, cat='Binary')

# Step 6: Define the objective function
# Minimize the total number of unicorns
prob += pulp.lpSum(unicorns[c] for c in cells), "Total_Unicorns"

# Step 7: Define the constraints
# For each black square, at least one unicorn must attack it.
for b_sq in black_squares:
    bx, by, bz = b_sq
    
    # The set of cells from which a unicorn can attack b_sq
    attacking_cells = []
    for s_sq in cells:
        sx, sy, sz = s_sq
        # A unicorn at s_sq attacks b_sq if the move is diagonal.
        # This includes placing a unicorn on the square itself (diffs are all 0).
        if abs(sx - bx) == abs(sy - by) == abs(sz - bz):
            attacking_cells.append(s_sq)
    
    # The sum of unicorns on attacking squares must be >= 1
    prob += pulp.lpSum(unicorns[s] for s in attacking_cells) >= 1, f"Cover_{b_sq}"

# Step 8: Solve the problem
prob.solve()

# Step 9: Print the results
print(f"Status: {pulp.LpStatus[prob.status]}")

if pulp.LpStatus[prob.status] == 'Optimal':
    min_unicorns = int(pulp.value(prob.objective))
    print(f"\nThe minimum number of unicorns needed is: {min_unicorns}")

    # As requested, output the final equation representing the sum.
    # For example: 1 + 1 + 1 + 1 = 4
    equation_parts = ["1"] * min_unicorns
    equation_str = " + ".join(equation_parts)
    print(f"The final equation is: {equation_str} = {min_unicorns}")

    # Find the positions of the placed unicorns for full detail
    # placed_unicorns = []
    # for c in cells:
    #     if unicorns[c].value() == 1.0:
    #         placed_unicorns.append(c)
    # print("\nOptimal unicorn placements:")
    # for pos in sorted(placed_unicorns):
    #     print(pos)
else:
    print("Could not find the optimal solution.")
