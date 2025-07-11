import pulp

def solve_unicorn_domination():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    # This script requires the 'pulp' library.
    # You can install it via pip: pip install pulp

    # 1. Define the board dimensions and identify cell types
    N = 8
    positions = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]
    # A cell (x, y, z) is black if the sum of its coordinates is odd.
    black_cells = [p for p in positions if (p[0] + p[1] + p[2]) % 2 != 0]

    # 2. Set up the ILP model
    # We want to minimize the number of unicorns, so it's an LpMinimize problem.
    model = pulp.LpProblem("Unicorn_Dominating_Set", pulp.LpMinimize)

    # 3. Define the variables
    # Create a binary variable for each cell on the board.
    # The variable is 1 if a unicorn is placed on the cell, 0 otherwise.
    unicorns = pulp.LpVariable.dicts("Unicorn", positions, cat='Binary')

    # 4. Define the objective function
    # The objective is to minimize the total sum of unicorns placed.
    model += pulp.lpSum([unicorns[p] for p in positions]), "Total_Unicorns"

    # 5. Define the constraints
    # Helper function to check if a unicorn at pos1 attacks pos2.
    # A unicorn moves along triagonals, where the change in x, y, and z is equal.
    def can_attack(pos1, pos2):
        x1, y1, z1 = pos1
        x2, y2, z2 = pos2
        # A unicorn on a square also attacks that square (distance = 0)
        return abs(x1 - x2) == abs(y1 - y2) == abs(z1 - z2)

    # For every black cell, we must ensure it's attacked by at least one unicorn.
    print("Setting up constraints for all 256 black squares... (This may take a moment)")
    for j in black_cells:
        # The sum of unicorns on all positions 'i' that can attack 'j' must be >= 1.
        model += pulp.lpSum([unicorns[i] for i in positions if can_attack(i, j)]) >= 1, f"Attack_Constraint_for_{j}"

    # 6. Solve the model
    print("Solving the Integer Linear Program...")
    # The solver will find the optimal solution. PuLP uses the CBC solver by default.
    model.solve()

    # 7. Print the results
    print("\n--- Solution ---")
    if pulp.LpStatus[model.status] == 'Optimal':
        min_unicorns = int(pulp.value(model.objective))
        
        print(f"The minimum number of unicorns required to attack all black squares is: {min_unicorns}")

        # Construct and print the equation as requested
        equation_parts = ["1"] * min_unicorns
        equation_str = " + ".join(equation_parts)
        print("\nThe equation for the minimum number is:")
        print(f"{equation_str} = {min_unicorns}")
        
        # Optionally, print the locations of the placed unicorns
        print("\nAn example of optimal unicorn placement is:")
        for p in positions:
            if unicorns[p].value() == 1:
                print(f"- {p}")
    else:
        print(f"Could not find the optimal solution. Status: {pulp.LpStatus[model.status]}")

# Execute the function
solve_unicorn_domination()