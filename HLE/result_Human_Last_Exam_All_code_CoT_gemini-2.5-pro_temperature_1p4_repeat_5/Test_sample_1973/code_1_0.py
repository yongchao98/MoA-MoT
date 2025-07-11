import pulp

def solve_unicorn_covering_problem():
    """
    This function formulates and solves the problem of finding the minimum
    number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    
    # Size of the 3D board
    N = 8

    # --- 1. Define the ILP Model ---
    # We want to minimize the number of unicorns
    prob = pulp.LpProblem("Unicorn_Covering_Black_Squares", pulp.LpMinimize)

    # --- 2. Define Board Coordinates and Decision Variables ---
    # List of all possible square coordinates (tuples)
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]
    
    # A square (x,y,z) is black if the sum of its coordinates is odd
    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 == 1]

    # Decision variables: Unicorn_(x,y,z) is 1 if a unicorn is on (x,y,z), 0 otherwise
    variables = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # --- 3. Define the Objective Function ---
    # Minimize the total number of placed unicorns
    prob += pulp.lpSum(variables[s] for s in all_squares), "Total_Unicorns"

    # --- 4. Define the Constraints ---
    # For each black square, we need a constraint to ensure it's attacked.
    for b_sq in black_squares:
        # The set of squares from which a unicorn can attack b_sq
        attacker_squares = []
        for s_sq in all_squares:
            # A unicorn move requires the absolute difference in each coordinate to be equal and non-zero.
            # We assume a piece does not attack the square it occupies.
            if s_sq == b_sq:
                continue
            
            dx = abs(s_sq[0] - b_sq[0])
            dy = abs(s_sq[1] - b_sq[1])
            dz = abs(s_sq[2] - b_sq[2])
            
            if dx > 0 and dx == dy and dx == dz:
                attacker_squares.append(s_sq)
        
        # Constraint: The sum of unicorns on all attacking squares must be at least 1.
        if attacker_squares:
            prob += pulp.lpSum(variables[s] for s in attacker_squares) >= 1, f"Attack_Constraint_for_{b_sq}"

    # --- 5. Solve the ILP Problem ---
    # We use the default CBC solver that comes with PuLP.
    # The `msg=False` argument suppresses solver output.
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # --- 6. Print the Solution ---
    if prob.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(prob.objective))
        print(f"The minimum number of unicorns to attack all black squares is: {min_unicorns}")
        
        # Find the locations of the placed unicorns
        unicorn_placements = [s for s in all_squares if variables[s].varValue > 0.5]
        
        # Construct the final equation string
        equation_parts = [f"1 * Unicorn{s}" for s in unicorn_placements]
        final_equation = " + ".join(equation_parts) + f" = {min_unicorns}"
        
        print("\nThe final equation representing the solution is:")
        print(final_equation)
        
        # Output the final answer in the requested format
        print("\n<<<{}>>>".format(min_unicorns))

    else:
        print("An optimal solution could not be found.")
        print(f"Solver status: {pulp.LpStatus[prob.status]}")

# Execute the function to find the solution
solve_unicorn_covering_problem()