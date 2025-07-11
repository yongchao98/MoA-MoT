import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 chessboard using Integer Linear Programming.
    """
    N = 8

    # Define all squares and identify the black ones
    all_squares = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    black_squares = [sq for sq in all_squares if sum(sq) % 2 != 0]

    # Create the ILP problem
    prob = pulp.LpProblem("Min_Unicorn_Cover", pulp.LpMinimize)

    # Decision Variables: A binary variable for each square, indicating if a unicorn is placed there
    unicorns = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # Objective Function: Minimize the total number of unicorns
    prob += pulp.lpSum(unicorns), "Total_Unicorns"

    # Constraints: Each black square must be attacked by at least one unicorn
    print("Generating constraints for each of the 256 black squares...")
    for b_sq in black_squares:
        r, c, l = b_sq
        
        # A square is attacked if a unicorn is on the same 3D diagonal
        # This includes the square itself (distance = 0)
        # The condition is abs(i-r) == abs(j-c) == abs(k-l)
        attacking_unicorns_sum = pulp.lpSum(
            unicorns[sq] for sq in all_squares
            if abs(sq[0] - r) == abs(sq[1] - c) and abs(sq[1] - c) == abs(sq[2] - l)
        )
        
        # Add the constraint for this black square
        prob += attacking_unicorns_sum >= 1, f"Attack_Coverage_for_{r}_{c}_{l}"

    # Solve the problem
    print("Solving the Integer Linear Program. This may take a moment...")
    # Using the default CBC solver from PuLP. It's effective for this scale.
    prob.solve()

    # --- Output the results ---
    print("\n" + "="*30)
    print("         SOLUTION")
    print("="*30)
    print(f"Solver Status: {pulp.LpStatus[prob.status]}")

    if prob.status == pulp.LpStatusOptimal:
        min_unicorns_count = int(pulp.value(prob.objective))
        print(f"\nThe minimum number of unicorns required is: {min_unicorns_count}")

        # Find where the unicorns are placed in the optimal solution
        placed_unicorns = [sq for sq in all_squares if pulp.value(unicorns[sq]) == 1]
        
        print("\nOne possible placement for the unicorns:")
        for pos in placed_unicorns:
            print(f"Unicorn at {pos}")
            
        # Display the final equation as requested
        print("\nFinal equation representing the count:")
        equation_str = " + ".join(["1"] * len(placed_unicorns))
        print(f"{equation_str} = {min_unicorns_count}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_unicorn_problem()
