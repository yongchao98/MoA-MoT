import pulp

def solve_unicorn_problem():
    """
    Solves the 3D chess unicorn covering problem using Integer Linear Programming.
    """
    # This script requires the 'pulp' library. You can install it using:
    # pip install pulp

    # --- Step 1: Define the subproblem ---
    # We solve for one of the four independent subproblems.
    # This subproblem consists of covering all (odd,odd,odd) squares, which are black,
    # using unicorns placed on either (odd,odd,odd) or (even,even,even) squares.
    
    print("Setting up the ILP for one of the four subproblems...")

    even_coords = [0, 2, 4, 6]
    odd_coords = [1, 3, 5, 7]
    
    # Potential locations for unicorns in this subproblem
    eee_squares = [(x, y, z) for x in even_coords for y in even_coords for z in even_coords]
    ooo_squares = [(x, y, z) for x in odd_coords for y in odd_coords for z in odd_coords]
    candidate_locations = eee_squares + ooo_squares

    # The black squares we need to cover in this subproblem
    target_squares = ooo_squares

    # --- Step 2: Create the ILP model ---
    model = pulp.LpProblem("Unicorn_Covering_Subproblem", pulp.LpMinimize)

    # --- Step 3: Define decision variables ---
    # A variable for each possible unicorn location
    unicorn_vars = pulp.LpVariable.dicts(
        "Unicorn", 
        candidate_locations, 
        cat='Binary'
    )

    # --- Step 4: Define the objective function ---
    # Minimize the total number of unicorns placed
    model += pulp.lpSum(unicorn_vars[loc] for loc in candidate_locations), "Total_Unicorns"

    # --- Step 5: Define the constraints ---
    # Each target black square must be attacked by at least one unicorn.
    print(f"Generating constraints for {len(target_squares)} black squares...")
    for target in target_squares:
        attackers = []
        for loc in candidate_locations:
            x1, y1, z1 = loc
            x2, y2, z2 = target
            # Unicorn attack rule: absolute differences in coords are equal
            if abs(x1 - x2) == abs(y1 - y2) == abs(z1 - z2):
                attackers.append(loc)
        
        model += pulp.lpSum(unicorn_vars[loc] for loc in attackers) >= 1, f"Cover_{target}"

    # --- Step 6: Solve the problem ---
    print("Solving the ILP problem... (This may take a moment)")
    # We use the default CBC solver that comes with PuLP, suppressing solver messages.
    model.solve(pulp.PULP_CBC_CMD(msg=0))

    # --- Step 7: Calculate and print the final result ---
    if model.status == pulp.LpStatusOptimal:
        min_unicorns_subproblem = int(pulp.value(model.objective))
        print(f"The minimum number of unicorns for one subproblem is: {min_unicorns_subproblem}")
        
        total_unicorns = 4 * min_unicorns_subproblem
        print("\nSince the board decomposes into 4 identical and independent subproblems, the total number is:")
        print(f"4 * {min_unicorns_subproblem} = {total_unicorns}")
    else:
        print("Could not find the optimal solution.")
        print(f"Status: {pulp.LpStatus[model.status]}")

if __name__ == '__main__':
    solve_unicorn_problem()
