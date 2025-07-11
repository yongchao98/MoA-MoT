import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 board using Integer Linear Programming.
    """
    N = 8

    # --- Step 1: Decompose the problem ---
    # The board can be split into 4 independent subproblems. We solve for one.
    # We choose the set where x, y, and z all have the same parity.
    # This corresponds to (x-y)%2 == 0 and (y-z)%2 == 0.
    
    subproblem_squares = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if (x % 2) == (y % 2) and (y % 2) == (z % 2):
                    subproblem_squares.append((x, y, z))

    # Identify the black squares to be covered in this subproblem.
    # A square (x,y,z) is black if (x+y+z) is odd.
    black_squares_to_cover = []
    for s in subproblem_squares:
        if sum(s) % 2 != 0:
            black_squares_to_cover.append(s)

    # --- Step 2: Formulate the ILP (Set Cover) ---
    
    # Create the model
    model = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # Create binary decision variables for placing a unicorn on any square in the set
    unicorn_vars = pulp.LpVariable.dicts(
        "Unicorn",
        subproblem_squares,
        cat='Binary'
    )

    # Define the objective function: minimize the number of unicorns
    model += pulp.lpSum(unicorn_vars[s] for s in subproblem_squares), "Total_Unicorns"

    # Define the constraints: each black square must be attacked
    for j_square in black_squares_to_cover:
        jx, jy, jz = j_square
        attackers = []
        for i_square in subproblem_squares:
            ix, iy, iz = i_square
            # A unicorn at i attacks j if the diagonal condition is met
            if abs(ix - jx) == abs(iy - jy) == abs(iz - jz):
                attackers.append(i_square)
        
        # The sum of unicorns on attacking squares must be at least 1
        model += pulp.lpSum(unicorn_vars[i] for i in attackers) >= 1, f"Cover_{j_square}"

    # --- Step 3: Solve the subproblem ---
    print("Solving the ILP for one of the four independent subproblems...")
    # Use the default CBC solver, suppress verbose output
    model.solve(pulp.PULP_CBC_CMD(msg=0))
    
    # --- Step 4: Output the results ---
    print(f"ILP solver status: {pulp.LpStatus[model.status]}")
    
    if model.status == pulp.LpStatusOptimal:
        min_unicorns_subproblem = int(pulp.value(model.objective))
        print(f"Minimum unicorns required for one subproblem: {min_unicorns_subproblem}")
        
        total_min_unicorns = 4 * min_unicorns_subproblem
        
        print("\nThe total minimum number of unicorns for the full 8x8x8 board is 4 times this value.")
        print("Final Equation:")
        print(f"{total_min_unicorns} = 4 * {min_unicorns_subproblem}")
        
        # The final answer in the required format
        print(f"\n<<<{total_min_unicorns}>>>")
    else:
        print("Could not find the optimal solution.")

solve_unicorn_problem()