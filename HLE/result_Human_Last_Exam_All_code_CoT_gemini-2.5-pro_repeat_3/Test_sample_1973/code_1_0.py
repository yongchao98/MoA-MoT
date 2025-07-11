import pulp

def solve_unicorn_dominance():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    N = 8
    
    # 1. Define all cells and identify the black cells
    print("Step 1: Defining the 8x8x8 board and identifying black squares...")
    all_cells = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]
    black_cells = [cell for cell in all_cells if sum(cell) % 2 != 0]
    print(f"Board has {len(all_cells)} cells in total.")
    print(f"There are {len(black_cells)} black squares to attack.\n")

    # 2. Determine which squares attack which black squares
    print("Step 2: Calculating all possible unicorn attacks on black squares...")
    attacks_map = {p: [] for p in all_cells}
    for unicorn_pos in all_cells:
        ux, uy, uz = unicorn_pos
        for black_cell in black_cells:
            bx, by, bz = black_cell
            # A unicorn attacks along space diagonals
            if abs(ux - bx) == abs(uy - by) == abs(uz - bz):
                attacks_map[unicorn_pos].append(black_cell)
    print("Attack map calculated.\n")

    # 3. Set up the ILP problem using PuLP
    print("Step 3: Formulating the Integer Linear Programming problem...")
    # Create the model
    prob = pulp.LpProblem("Unicorn_Cover", pulp.LpMinimize)

    # Decision Variables: A binary variable for each cell, 1 if a unicorn is placed, 0 otherwise
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", all_cells, cat=pulp.LpBinary)

    # Objective Function: Minimize the total number of unicorns
    prob += pulp.lpSum([unicorn_vars[p] for p in all_cells]), "Total_Unicorns"

    # Constraints: Each black square must be covered by at least one unicorn
    for b_cell in black_cells:
        # Find all possible unicorn positions that attack this black cell
        possible_attackers = [p for p in all_cells if b_cell in attacks_map[p]]
        # Add the constraint
        prob += pulp.lpSum([unicorn_vars[p] for p in possible_attackers]) >= 1, f"Cover_{b_cell}"
    print("ILP model created.\n")

    # 4. Solve the ILP
    print("Step 4: Solving the ILP problem (this may take a moment)...")
    prob.solve()
    print("Solver finished.\n")

    # 5. Print the results
    print("--- Solution ---")
    print(f"Status: {pulp.LpStatus[prob.status]}")
    
    min_unicorns = int(pulp.value(prob.objective))
    print(f"Minimum number of unicorns required: {min_unicorns}")

    placed_unicorns = [p for p in all_cells if unicorn_vars[p].value() == 1.0]

    # Format the final equation as requested
    equation_parts = [f"1*U{p}" for p in placed_unicorns]
    print("\nFinal Equation:")
    print(f"{' + '.join(equation_parts)} = {min_unicorns}")

    print("\nWhere the unicorn positions U(x,y,z) are:")
    for i, p in enumerate(placed_unicorns):
        print(f"  Unicorn {i+1}: {p}")

if __name__ == '__main__':
    solve_unicorn_dominance()
    print("\n<<<32>>>")