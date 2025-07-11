import pulp

def solve_subproblem(problem_name, target_parity, attacker_parities):
    """
    Solves a subproblem for a specific group of squares defined by coordinate parity.

    Args:
        problem_name (str): The name for the ILP problem.
        target_parity (tuple): A 3-tuple of parities (0 for even, 1 for odd) for the squares to be attacked.
        attacker_parities (list): A list of parity tuples for squares where unicorns can be placed.

    Returns:
        int: The minimum number of unicorns required for this subproblem.
    """
    N = 8
    prob = pulp.LpProblem(problem_name, pulp.LpMinimize)

    # Define all coordinates and categorize them by their parity type
    all_coords = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    coord_parities = {c: (c[0] % 2, c[1] % 2, c[2] % 2) for c in all_coords}

    # Filter coordinates for target squares and potential unicorn locations
    target_squares = [c for c, p in coord_parities.items() if p == target_parity]
    attacker_candidate_squares = [c for c, p in coord_parities.items() if p in attacker_parities]

    # Create binary decision variables for each potential unicorn placement
    unicorns = pulp.LpVariable.dicts("unicorn", attacker_candidate_squares, cat=pulp.LpBinary)

    # Set the objective function to minimize the number of placed unicorns
    prob += pulp.lpSum(unicorns[c] for c in attacker_candidate_squares), "Total_Unicorns"

    # Add constraints: each target square must be attacked by at least one unicorn
    for b_coord in target_squares:
        bx, by, bz = b_coord
        # Sum of unicorns that attack the current target square
        attacking_sum = pulp.lpSum(
            unicorns[a_coord] for a_coord in attacker_candidate_squares
            if abs(bx - a_coord[0]) == abs(by - a_coord[1]) == abs(bz - a_coord[2])
        )
        prob += attacking_sum >= 1, f"Constraint_Attack_{bx}_{by}_{bz}"
    
    # Solve the ILP problem quietly
    prob.solve(pulp.PULP_CBC_CMD(msg=False))
    
    # Return the optimal value found
    return int(pulp.value(prob.objective))

def solve_raumschach_unicorn_cover():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board.
    """
    print("Solving subproblems. This may take a moment...")

    # Define the parity groups for black squares and their corresponding attacker groups
    # Using 0 for even, 1 for odd. A square (x,y,z) is black if x+y+z is odd (with 0-based indexing).
    
    # Subproblem Type 1: Symmetrical groups
    # e.g., covering black squares of type (even, even, odd)
    p_black_symm = (0, 0, 1) 
    # These can be attacked by unicorns on (e,e,o) or (o,o,e) squares
    p_attack_symm = [(0, 0, 1), (1, 1, 0)]
    n1 = solve_subproblem("Symmetric_Group_Cover", p_black_symm, p_attack_symm)
    print(f"Minimum unicorns for one of the three symmetric black square groups: {n1}")

    # Subproblem Type 2: The unique group
    # Covering black squares of type (odd, odd, odd)
    p_black_unique = (1, 1, 1)
    # These can be attacked by unicorns on (o,o,o) or (e,e,e) squares
    p_attack_unique = [(1, 1, 1), (0, 0, 0)]
    n2 = solve_subproblem("Unique_Group_Cover", p_black_unique, p_attack_unique)
    print(f"Minimum unicorns for the (odd,odd,odd) black square group: {n2}")
    
    # Calculate and print the total result
    total_unicorns = 3 * n1 + n2
    
    print("\n--- Final Result ---")
    print(f"The minimum number of unicorns needed is the sum of the solutions to the four subproblems.")
    print("The final equation is:")
    
    equation_str = f"{n1} + {n1} + {n1} + {n2} = {total_unicorns}"
    print(equation_str)
    
    print(f"\nThe minimum number of unicorns to attack all black squares is: {total_unicorns}")
    return total_unicorns

if __name__ == '__main__':
    min_unicorns_result = solve_raumschach_unicorn_cover()
    # The final line below is for the answer format extraction.
    print(f"\n<<<{min_unicorns_result}>>>")
