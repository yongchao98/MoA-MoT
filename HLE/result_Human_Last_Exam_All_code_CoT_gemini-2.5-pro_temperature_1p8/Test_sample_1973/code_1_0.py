import pulp

def solve_unicorn_attack_problem():
    """
    This function formulates and solves the problem of finding the minimum number
    of unicorns to attack all black squares on an 8x8x8 board using Integer
    Linear Programming.
    """
    N = 8

    # Step 1: Define the board squares and identify black squares.
    # Coordinates are 0-indexed, from 0 to 7.
    squares = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    # A square is black if the sum of its coordinates is odd.
    black_squares = [s for s in squares if (s[0] + s[1] + s[2]) % 2 != 0]

    # Step 2: Formulate the ILP problem using pulp.
    # We want to minimize the number of unicorns.
    model = pulp.LpProblem("Unicorn_Cover_Problem", pulp.LpMinimize)

    # Decision Variables: A binary variable for each square.
    # unicorn_vars[s] is 1 if a unicorn is at square s, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", squares, cat='Binary')

    # Objective Function: Minimize the total number of unicorns.
    model += pulp.lpSum(unicorn_vars[s] for s in squares)

    # Step 3: Define the constraints. Every black square must be attacked.
    # For efficiency, pre-calculate the set of squares that attack each square.
    print("Pre-calculating attack patterns...")
    attackers_map = {s: set() for s in squares}
    for s1 in squares:
        x1, y1, z1 = s1
        for s2 in squares:
            x2, y2, z2 = s2
            # A unicorn move means absolute differences in coordinates are equal.
            dx = abs(x1 - x2)
            dy = abs(y1 - y2)
            dz = abs(z1 - z2)
            if dx == dy and dy == dz:
                # If s2 can attack s1, add it to s1's list of attackers.
                attackers_map[s1].add(s2)
    print("Calculation complete.")

    # Add a constraint for each black square.
    for b_square in black_squares:
        # The sum of unicorns on attacking squares must be at least 1.
        model += pulp.lpSum(unicorn_vars[s] for s in attackers_map[b_square]) >= 1, f"Attack_Constraint_for_{b_square}"

    # Step 4: Solve the ILP problem.
    # This may take a few minutes as the problem is complex.
    print("Solving the ILP problem... This might take a few minutes.")
    model.solve()
    print("Solver finished.")

    # Step 5: Output the result.
    print("\n--- Results ---")
    status = pulp.LpStatus[model.status]
    print(f"Solution Status: {status}")

    if status == 'Optimal':
        min_unicorns = int(pulp.value(model.objective))
        print(f"The minimum number of unicorns required is: {min_unicorns}")

        # Find the locations of the unicorns (optional, for verification)
        # unicorn_placements = [s for s in squares if unicorn_vars[s].varValue == 1]
        # print(f"Found {len(unicorn_placements)} unicorn placements.")
        
        # As requested, output the final equation.
        equation_parts = ["1"] * min_unicorns
        print("\nFinal Equation:")
        print(f"{' + '.join(equation_parts)} = {min_unicorns}")
    else:
        print("An optimal solution could not be found.")

if __name__ == '__main__':
    solve_unicorn_attack_problem()