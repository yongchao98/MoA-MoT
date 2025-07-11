import pulp

def solve_unicorn_domination():
    """
    Solves the 3D chess unicorn problem using Integer Linear Programming.
    Finds the minimum number of unicorns to attack all black squares on an 8x8x8 board.
    """
    N = 8
    
    # This problem may take a minute or two to solve depending on your system.
    print(f"Setting up the model for an {N}x{N}x{N} board...")

    # 1. Define all squares and identify the black ones.
    # We use coordinates from 1 to N to match the coloring rule (x+y+z is even).
    all_squares = []
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            for k in range(1, N + 1):
                all_squares.append((i, j, k))

    black_squares = []
    for square in all_squares:
        if (square[0] + square[1] + square[2]) % 2 == 0:
            black_squares.append(square)
            
    print(f"Total squares: {len(all_squares)}")
    print(f"Black squares to cover: {len(black_squares)}")

    # 2. Create the ILP model
    model = pulp.LpProblem("Unicorn_Domination", pulp.LpMinimize)

    # 3. Define the decision variables
    # u_ijk = 1 if a unicorn is placed on square (i,j,k), 0 otherwise
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # 4. Set the objective function: Minimize the total number of unicorns
    model += pulp.lpSum(unicorn_vars[s] for s in all_squares), "Total_Unicorns"

    # 5. Add the constraints
    # For each black square, at least one unicorn must attack it.
    print("Generating constraints...")
    for b_square in black_squares:
        r, s, t = b_square
        
        # Find all squares p that can attack b_square
        # A square p=(i,j,k) attacks b=(r,s,t) if they are on the same 3D diagonal.
        # This means abs(i-r) == abs(j-s) == abs(k-t)
        attacking_positions = []
        for p_square in all_squares:
            i, j, k = p_square
            if abs(i - r) == abs(j - s) and abs(j - s) == abs(k - t):
                attacking_positions.append(p_square)
        
        # Add the constraint for this black square
        model += pulp.lpSum(unicorn_vars[p] for p in attacking_positions) >= 1, f"Attack_Coverage_for_{b_square}"

    # 6. Solve the model
    print("Solving the ILP problem... (this might take a moment)")
    # Using the default CBC solver. If you have Gurobi or CPLEX, it might be faster.
    model.solve()

    # 7. Print the results
    status = pulp.LpStatus[model.status]
    min_unicorns = int(pulp.value(model.objective))

    print("\n--- Solution ---")
    print(f"Status: {status}")
    
    if status == 'Optimal':
        print(f"The minimum number of unicorns required is: {min_unicorns}")
        
        # Print the equation as requested
        print("\nFinal equation:")
        equation = " + ".join(["1"] * min_unicorns)
        print(f"{equation} = {min_unicorns}")
        
        # This is the final numerical answer for the platform.
        print(f"\n<<<answer content>>>\n{min_unicorns}")
    else:
        print("Could not find the optimal solution.")


if __name__ == '__main__':
    # You might need to install pulp: pip install pulp
    solve_unicorn_domination()
