import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 board using Integer Linear Programming.
    """
    N = 8

    # 1. Generate and classify all squares on the board
    white_squares = []
    black_squares = []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                square = (i, j, k)
                if (i + j + k) % 2 == 0:
                    white_squares.append(square)
                else:
                    black_squares.append(square)

    # 2. Create the ILP problem
    prob = pulp.LpProblem("Unicorn_Dominating_Set", pulp.LpMinimize)

    # 3. Create decision variables: one for each white square
    # The keys are the square coordinates (tuples)
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", white_squares, cat='Binary')

    # 4. Set the objective function: Minimize the total number of unicorns
    prob += pulp.lpSum(unicorn_vars[w] for w in white_squares), "Total_Unicorns"

    # 5. Add constraints: Each black square must be attacked by at least one unicorn
    for b_square in black_squares:
        b_i, b_j, b_k = b_square
        
        # Find all white squares that can attack the current black square
        attacking_unicorns = []
        for w_square in white_squares:
            w_i, w_j, w_k = w_square
            
            # Check the unicorn's attack condition
            di = abs(b_i - w_i)
            dj = abs(b_j - w_j)
            dk = abs(b_k - w_k)
            
            if di > 0 and di == dj and di == dk:
                attacking_unicorns.append(unicorn_vars[w_square])
        
        # Add the constraint for this black square. The sum of unicorns
        # on attacking squares must be at least 1.
        if attacking_unicorns:
            prob += pulp.lpSum(attacking_unicorns) >= 1, f"Attack_on_{b_i}_{b_j}_{b_k}"

    # 6. Solve the problem (suppress solver messages)
    # This might take a minute or two depending on your system's performance.
    print("Solving the Integer Linear Programming problem... (This may take a moment)")
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    print("Solver finished.")

    # 7. Print the results
    if prob.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(prob.objective))
        print("\n--- Solution ---")
        print(f"The minimum number of unicorns needed to attack all black squares is: {min_unicorns}")
        
        # Create and print the equation as requested
        print("\nThe equation representing the sum of unicorns is:")
        equation_parts = ["1"] * min_unicorns
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {min_unicorns}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_unicorn_problem()
