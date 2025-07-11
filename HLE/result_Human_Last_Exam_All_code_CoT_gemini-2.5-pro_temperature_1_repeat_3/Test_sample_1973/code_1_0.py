import pulp

def solve_unicorn_domination():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    N = 8  # Board size

    # Generate all possible squares (cells) on the 8x8x8 board
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

    # Identify the black squares. A square (x, y, z) is black if x+y+z is odd.
    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 != 0]

    # --- ILP Model Setup ---

    # 1. Create the minimization problem
    prob = pulp.LpProblem("Unicorn_Black_Square_Cover", pulp.LpMinimize)

    # 2. Define Decision Variables
    # u_vars[square] = 1 if a unicorn is placed at 'square', 0 otherwise.
    u_vars = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # 3. Define the Objective Function
    # Minimize the total number of placed unicorns.
    prob += pulp.lpSum(u_vars[s] for s in all_squares), "Total_Unicorns"

    # 4. Define the Constraints
    # For each black square, at least one unicorn must attack it.
    for b_sq in black_squares:
        i, j, k = b_sq
        
        # A unicorn at (x,y,z) attacks (i,j,k) if they are on a 3D diagonal.
        # This is true if |x-i| = |y-j| = |z-k|.
        # This is equivalent to checking four families of invariants.
        attacking_unicorns = []
        for p_sq in all_squares:
            x, y, z = p_sq
            # Check the four diagonal conditions
            cond1 = ((y - x) == (j - i) and (z - x) == (k - i))
            cond2 = ((y - x) == (j - i) and (z + x) == (k + i))
            cond3 = ((y + x) == (j + i) and (z - x) == (k - i))
            cond4 = ((y + x) == (j + i) and (z + x) == (k + i))

            if cond1 or cond2 or cond3 or cond4:
                attacking_unicorns.append(u_vars[p_sq])
        
        # Add the constraint: the sum of attacking unicorns must be >= 1
        prob += pulp.lpSum(attacking_unicorns) >= 1, f"Constraint_for_square_{i}_{j}_{k}"

    # --- Solve the problem ---
    prob.solve()

    # --- Print the results ---
    print(f"Solver status: {pulp.LpStatus[prob.status]}")
    
    min_unicorns = int(pulp.value(prob.objective))
    print(f"The minimum number of unicorns to attack all black squares is: {min_unicorns}")
    
    print("\nAn example of an optimal placement is:")
    unicorn_positions = [s for s in all_squares if u_vars[s].varValue > 0.5]
    for pos in unicorn_positions:
        print(f"Unicorn at {pos}")

    # As requested, output the final equation representing the sum.
    equation_str = " + ".join(["1"] * min_unicorns)
    print(f"\nThe final equation is: {equation_str} = {min_unicorns}")

if __name__ == '__main__':
    solve_unicorn_domination()