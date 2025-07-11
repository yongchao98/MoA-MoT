import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board.
    """
    N = 8
    
    # --- 1. Define the Board and Squares ---
    
    # Create a list of all cell coordinates (from 0 to 7)
    cells = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    
    # Identify the black squares. A square (i, j, k) is black if (i+j+k) is even.
    black_cells = [c for c in cells if (c[0] + c[1] + c[2]) % 2 == 0]

    # --- 2. Determine Attack Moves ---

    # Pre-compute which cells can attack which other cells to speed up constraint creation.
    # A unicorn at c1 attacks c2 if they are on a space diagonal.
    # This means the absolute difference in each coordinate is the same.
    # A unicorn also "attacks" the square it sits on (distance = 0).
    attack_map = {target: [] for target in cells}
    for c1 in cells:
        x1, y1, z1 = c1
        for c2 in cells:
            x2, y2, z2 = c2
            if abs(x1 - x2) == abs(y1 - y2) == abs(z1 - z2):
                attack_map[c2].append(c1)

    # --- 3. Formulate the ILP Problem ---
    
    # We want to MINIMIZE the number of unicorns.
    model = pulp.LpProblem("Unicorn_Dominating_Set", pulp.LpMinimize)

    # Define a binary variable for each cell on the board.
    # The variable is 1 if a unicorn is placed there, 0 otherwise.
    unicorns = pulp.LpVariable.dicts("Unicorn_Placement", cells, cat=pulp.LpBinary)

    # The objective is to minimize the total number of placed unicorns.
    model += pulp.lpSum([unicorns[c] for c in cells]), "Total_Unicorns"

    # Add a constraint for each black cell: it must be attacked by at least one unicorn.
    for b_cell in black_cells:
        # The sum of unicorns on squares that can attack this black cell must be >= 1.
        model += pulp.lpSum([unicorns[attacker] for attacker in attack_map[b_cell]]) >= 1, f"Constraint_for_{b_cell}"

    # --- 4. Solve the Problem ---
    
    # Solve the model. Pulp will call a solver (like CBC) to find the optimal solution.
    # This may take a moment.
    print("Solving the Integer Linear Program... This might take a few minutes.")
    model.solve()
    print("Solver finished.")

    # --- 5. Display the Result ---
    
    if model.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(model.objective))
        print("\n--- Solution ---")
        print(f"Status: Optimal")
        print(f"Minimum number of unicorns to attack all black squares: {min_unicorns}")
        
        # Fulfilling the request to show the final equation
        equation_parts = ["1"] * min_unicorns
        equation = " + ".join(equation_parts)
        if min_unicorns > 0:
            print(f"Equation: {equation} = {min_unicorns}")
        else:
            print("Equation: 0 = 0")
        
        # Returning the final answer in the specified format
        print(f"\n<<<__{min_unicorns}__>>>")

    else:
        print(f"Could not find the optimal solution. Status: {pulp.LpStatus[model.status]}")

if __name__ == '__main__':
    solve_unicorn_problem()