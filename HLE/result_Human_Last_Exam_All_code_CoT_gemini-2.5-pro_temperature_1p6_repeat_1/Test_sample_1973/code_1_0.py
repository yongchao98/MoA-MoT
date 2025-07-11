import pulp

def solve_unicorn_covering_problem():
    """
    Finds the minimum number of unicorns to attack all black squares on an 8x8x8
    chessboard using Integer Linear Programming.
    """
    BOARD_SIZE = 8

    # 1. Define the ILP model. It's a minimization problem.
    model = pulp.LpProblem("Unicorn_Covering_Problem", pulp.LpMinimize)

    # 2. Define all squares and identify the black ones.
    all_squares = [(x, y, z) for x in range(BOARD_SIZE) for y in range(BOARD_SIZE) for z in range(BOARD_SIZE)]
    black_squares = [s for s in all_squares if sum(s) % 2 != 0]

    # 3. Create a binary decision variable for each square.
    # U[s] = 1 if a unicorn is placed on square s, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts("Unicorn_at", all_squares, cat='Binary')

    # 4. Set the objective function: Minimize the total number of unicorns.
    model += pulp.lpSum(unicorn_vars[s] for s in all_squares), "Minimize_Total_Unicorns"

    # 5. Add constraints: Each black square must be attacked by at least one unicorn.
    for b_sq in black_squares:
        bx, by, bz = b_sq
        
        # Find all squares s_sq from which a unicorn can attack the black square b_sq.
        attacking_squares = []
        for s_sq in all_squares:
            sx, sy, sz = s_sq
            # The attack condition for a unicorn move.
            if abs(sx - bx) == abs(sy - by) and abs(sy - by) == abs(sz - bz):
                attacking_squares.append(s_sq)
        
        # Add the constraint that the sum of unicorns on attacking squares must be >= 1.
        model += pulp.lpSum(unicorn_vars[s] for s in attacking_squares) >= 1, f"Cover_Black_Square_{b_sq}"

    # 6. Solve the problem. This might take a few moments.
    print("Solving the ILP model to find the minimum number of unicorns...")
    model.solve()
    print("Solver has finished.")

    # 7. Print the results.
    if model.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(model.objective))
        print(f"\nMinimum number of unicorns needed to attack all black squares is: {min_unicorns}")

        placed_unicorns = [s for s in all_squares if unicorn_vars[s].varValue == 1]
        
        print("\nAn optimal placement solution is given by the following equation:")
        print("Let U_(x,y,z) be 1 if a unicorn is at (x,y,z), and 0 otherwise.")
        
        # Print each term of the final equation.
        for i, pos in enumerate(placed_unicorns):
            # The name of each term, e.g., U_(0,1,2)
            term = f"U_{pos}"
            print(term, end="")
            if i < len(placed_unicorns) - 1:
                print(" + ", end="")
        
        print(f" = {min_unicorns}")
    else:
        print("\nCould not find the optimal solution.")
        print(f"Solver Status: {pulp.LpStatus[model.status]}")

if __name__ == '__main__':
    solve_unicorn_covering_problem()