import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 board using Integer Linear Programming.
    """
    # This script requires the 'pulp' library.
    # You can install it with: pip install pulp
    
    # 1. Define the Board and Squares
    N = 8
    all_squares = []
    black_squares = []

    # A square (x,y,z) is black if (x+y+z) is odd.
    print("Step 1: Identifying all squares and black squares on the 8x8x8 board...")
    for x in range(N):
        for y in range(N):
            for z in range(N):
                pos = (x, y, z)
                all_squares.append(pos)
                if (x + y + z) % 2 == 1:
                    black_squares.append(pos)
    
    print(f"Total squares: {len(all_squares)}")
    print(f"Black squares to cover: {len(black_squares)}\n")

    # 2. Define the Unicorn's Attack Logic
    def can_attack(pos1, pos2):
        """Checks if a unicorn at pos1 can attack pos2."""
        x1, y1, z1 = pos1
        x2, y2, z2 = pos2
        # The move is along a space diagonal, so coordinate differences must be equal.
        # A unicorn also attacks its own square (distance k=0).
        return abs(x1 - x2) == abs(y1 - y2) == abs(z1 - z2)

    # 3. Formulate the ILP Problem
    print("Step 2: Setting up the Integer Linear Programming model...")
    
    # Create the minimization problem
    prob = pulp.LpProblem("Min_Unicorn_Cover", pulp.LpMinimize)

    # Create a binary variable for each square on the board.
    # x_pos = 1 if a unicorn is placed at pos, 0 otherwise.
    variables = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # Define the objective function: Minimize the total number of unicorns.
    prob += pulp.lpSum(variables[pos] for pos in all_squares), "TotalUnicorns"

    # Define the constraints: Each black square must be attacked.
    print("Step 3: Adding constraints - ensuring every black square is attacked...")
    for target_pos in black_squares:
        # Find all squares from which a unicorn can attack the target_pos
        attackers = [attacker_pos for attacker_pos in all_squares if can_attack(attacker_pos, target_pos)]
        # The sum of variables for squares that can attack the target must be >= 1.
        prob += pulp.lpSum(variables[pos] for pos in attackers) >= 1, f"Coverage_for_{target_pos}"

    # 4. Solve the Problem
    print("Step 4: Solving the ILP problem (this may take a moment)...")
    # CBC is the default solver that comes with PuLP
    prob.solve()
    print("\n" + "="*30)
    print("            SOLUTION")
    print("="*30)

    # 5. Output the Results
    status = pulp.LpStatus[prob.status]
    min_unicorns = int(pulp.value(prob.objective))

    print(f"Status: {status}")
    
    if status == 'Optimal':
        # Print the final equation as requested
        print("\nThe final equation for the minimum number is:")
        equation_parts = ["1"] * min_unicorns
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {min_unicorns}\n")
        
        print(f"The minimum number of unicorns needed to attack all black squares is: {min_unicorns}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_unicorn_problem()