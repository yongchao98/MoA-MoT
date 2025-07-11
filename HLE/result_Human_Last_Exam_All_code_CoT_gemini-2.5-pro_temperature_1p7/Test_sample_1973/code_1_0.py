import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8
    chessboard using Integer Linear Programming.
    """
    # This function requires the 'pulp' library.
    # You can install it with: pip install pulp

    # 1. Define board size and generate squares
    N = 8
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

    # 2. Classify squares into white and black
    # White squares: x+y+z is even. Unicorns will be placed here.
    # Black squares: x+y+z is odd. These are the squares to be attacked.
    white_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 == 0]
    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 == 1]

    # 3. Define the ILP problem
    prob = pulp.LpProblem("Unicorn_Covering_Problem", pulp.LpMinimize)

    # 4. Create decision variables
    # u_w = 1 if a unicorn is placed on white square w, 0 otherwise.
    # PuLP variables are stored in a dictionary. Keys are strings representing the coordinates.
    unicorn_vars = pulp.LpVariable.dicts(
        "Unicorn",
        [str(s) for s in white_squares],
        cat='Binary'
    )

    # 5. Set the objective function
    # Minimize the total number of placed unicorns.
    prob += pulp.lpSum(unicorn_vars[str(s)] for s in white_squares), "Total_Unicorns"

    # 6. Add constraints
    # Each black square must be attacked by at least one unicorn.
    for b_square in black_squares:
        bx, by, bz = b_square
        
        # Find all white squares that can attack this black square
        attacking_white_squares = []
        for w_square in white_squares:
            wx, wy, wz = w_square
            # Unicorn move condition: absolute differences in coordinates are equal
            if abs(wx - bx) == abs(wy - by) == abs(wz - bz):
                attacking_white_squares.append(w_square)
        
        # The sum of unicorns on attacking squares must be >= 1
        prob += (
            pulp.lpSum(unicorn_vars[str(s)] for s in attacking_white_squares) >= 1,
            f"Attack_Constraint_for_{b_square}"
        )

    # 7. Solve the ILP problem
    # The solver is automatically called by PuLP. By default, it uses CBC.
    # You can suppress the solver's own logs by using `prob.solve(pulp.PULP_CBC_CMD(msg=0))`
    prob.solve()

    # 8. Print the results
    print(f"Solver status: {pulp.LpStatus[prob.status]}")
    
    # Extract the final number of unicorns
    min_unicorns = int(pulp.value(prob.objective))
    print(f"The minimum number of unicorns needed is: {min_unicorns}")
    
    # As requested, output the final sum equation.
    # This will be a sum of '1' for each unicorn in the minimal solution.
    equation_str = " + ".join(["1"] * min_unicorns)
    print("\nFinal equation:")
    print(f"{equation_str} = {min_unicorns}")

if __name__ == '__main__':
    solve_unicorn_problem()

<<<16>>>