import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 3D chessboard using Integer Linear Programming.
    """
    # Note: This script requires the 'pulp' library.
    # You can install it using: pip install pulp

    board_size = 8

    # 1. Identify all white and black squares
    white_squares = []
    black_squares = []
    for x in range(board_size):
        for y in range(board_size):
            for z in range(board_size):
                if (x + y + z) % 2 == 0:
                    white_squares.append((x, y, z))
                else:
                    black_squares.append((x, y, z))

    # 2. Define the ILP problem
    prob = pulp.LpProblem("Unicorn_Dominating_Set", pulp.LpMinimize)

    # 3. Create a binary decision variable for each white square
    # This variable is 1 if we place a unicorn on that square, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts(
        "Unicorn_Placement",
        white_squares,
        cat='Binary'
    )

    # 4. Set the objective function: Minimize the total number of unicorns
    prob += pulp.lpSum(unicorn_vars[w] for w in white_squares), "Total_Unicorns"

    # 5. Add constraints: Each black square must be attacked by at least one unicorn
    for b_sq in black_squares:
        xb, yb, zb = b_sq
        
        # Find all white squares that can attack this black square
        attackers = []
        for w_sq in white_squares:
            xw, yw, zw = w_sq
            dx = abs(xb - xw)
            dy = abs(yb - yw)
            dz = abs(zb - zw)
            
            if dx > 0 and dx == dy and dx == dz:
                attackers.append(w_sq)
        
        # The sum of unicorn variables on attacking squares must be >= 1
        if attackers:
            prob += pulp.lpSum(unicorn_vars[w] for w in attackers) >= 1, f"Cover_Square_{xb}_{yb}_{zb}"

    # 6. Solve the ILP problem
    # The default CBC solver bundled with pulp is used.
    prob.solve()

    # 7. Print the results
    min_unicorns = int(pulp.value(prob.objective))
    
    print(f"Status: {pulp.LpStatus[prob.status]}")
    print(f"Minimum number of unicorns needed: {min_unicorns}")
    
    # As requested, output the final equation
    # This will be a sum of 1s, representing the unicorns placed.
    equation_parts = ['1'] * min_unicorns
    equation_str = " + ".join(equation_parts)
    print("The final equation for the total number of unicorns is:")
    print(f"{equation_str} = {min_unicorns}")

solve_unicorn_problem()
<<<16>>>