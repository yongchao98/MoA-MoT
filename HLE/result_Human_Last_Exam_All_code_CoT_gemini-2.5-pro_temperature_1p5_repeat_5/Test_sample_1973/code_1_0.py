import pulp

def solve_unicorn_covering():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    N = 8  # Board size

    # Generate all cells and identify the black ones
    cells = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    black_cells = [c for c in cells if (c[0] + c[1] + c[2]) % 2 == 1]

    # --- ILP Setup using PuLP ---

    # 1. Initialize the problem
    # We want to minimize the number of unicorns.
    prob = pulp.LpProblem("Unicorn_Covering_Problem", pulp.LpMinimize)

    # 2. Define Decision Variables
    # U_i_j_k = 1 if a unicorn is on cell (i,j,k), 0 otherwise
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", cells, cat='Binary')

    # 3. Define the Objective Function
    # Minimize the sum of all unicorn variables.
    prob += pulp.lpSum([unicorn_vars[c] for c in cells]), "Total_Unicorns"

    # 4. Define the Constraints
    # Each black cell must be attacked by at least one unicorn.
    print("Setting up constraints... (This may take a moment)")
    for b_cell in black_cells:
        r, c, l = b_cell
        # Find all cells from which a unicorn can attack the current black cell
        possible_attackers = []
        for a_cell in cells:
            i, j, k = a_cell
            if i == r and j == c and l == k:
                continue  # A unicorn doesn't attack its own cell

            dr = abs(r - i)
            dc = abs(c - j)
            dl = abs(l - k)

            if dr > 0 and dr == dc and dr == dl:
                possible_attackers.append(unicorn_vars[a_cell])
        
        # The constraint: sum of unicorns that can attack this cell must be >= 1
        if possible_attackers:
            prob += pulp.lpSum(possible_attackers) >= 1, f"Cover_({r},{c},{l})"

    # 5. Solve the Problem
    # The default CBC solver in PuLP is sufficient here.
    # msg=False suppresses solver output.
    print("Solving... (This can take a minute or two depending on your system)")
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # --- Output the Results ---
    print("\n--- Solution ---")
    status = pulp.LpStatus[prob.status]
    print(f"Status: {status}")

    if status == 'Optimal':
        min_unicorns = int(pulp.value(prob.objective))
        print(f"\nThe minimum number of unicorns needed to attack all black squares is: {min_unicorns}")

        # Find the coordinates of the placed unicorns
        unicorn_locations = []
        for c in cells:
            if pulp.value(unicorn_vars[c]) == 1:
                unicorn_locations.append(c)

        # Print the final "equation" as requested
        print("\nAn equation representing one possible optimal placement is:")
        
        # Sorting locations makes the output consistent
        unicorn_locations.sort()
        
        equation_parts = [f"U{loc}" for loc in unicorn_locations]
        print(" + ".join(equation_parts) + f" = {min_unicorns}")
        
        print("\nWhere U(x,y,z) represents a unicorn placed at coordinate (x,y,z). The coordinates are:")
        for loc in unicorn_locations:
            print(loc)
    else:
        print("Could not find an optimal solution.")
        
    return pulp.value(prob.objective)


if __name__ == '__main__':
    final_answer = solve_unicorn_covering()
    if final_answer is not None:
        print(f"\n<<<__{int(final_answer)}__>>>")
