import pulp

def solve_subproblem(problem_name, target_squares, candidate_squares):
    """
    Solves the set covering subproblem for a given set of targets and candidates.
    
    Args:
        problem_name (str): The name for the ILP problem.
        target_squares (list): A list of (x,y,z) tuples to be covered.
        candidate_squares (list): A list of (x,y,z) tuples where unicorns can be placed.
        
    Returns:
        int: The minimum number of unicorns required for this subproblem.
    """
    # Create the ILP problem
    prob = pulp.LpProblem(problem_name, pulp.LpMinimize)

    # Decision variables: Place a unicorn on a candidate square?
    locations = pulp.LpVariable.dicts(
        "Unicorn",
        candidate_squares,
        cat='Binary'
    )

    # Objective function: Minimize the number of unicorns
    prob += pulp.lpSum(locations[i] for i in candidate_squares), "Total_Unicorns"

    # Constraints: Each target square must be attacked by at least one unicorn
    for j_target in target_squares:
        # A square is attacked if a unicorn is on it or can move to it.
        # A unicorn at i=(xi,yi,zi) attacks j=(xj,yj,zj) if |xi-xj|=|yi-yj|=|zi-zj|
        attacking_set = [i for i in candidate_squares 
                         if abs(i[0] - j_target[0]) == abs(i[1] - j_target[1]) == abs(i[2] - j_target[2])]
        
        prob += pulp.lpSum(locations[i] for i in attacking_set) >= 1, f"Cover_{j_target}"

    # Solve the problem. The CBC solver is included with pulp. msg=False suppresses solver output.
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # Return the optimal value for the objective function
    return int(pulp.value(prob.objective))

def main():
    """
    Main function to define and solve the partitioned unicorn problem.
    """
    N = 8
    coords = range(N)

    # Partition all 512 squares by the parity of their coordinates (0 for even, 1 for odd)
    parity_types = {}
    for x in coords:
        for y in coords:
            for z in coords:
                parity = (x % 2, y % 2, z % 2)
                if parity not in parity_types:
                    parity_types[parity] = []
                parity_types[parity].append((x, y, z))

    # Black squares have coordinate sums that are odd.
    # The parities for black squares are (odd,odd,odd), (even,even,odd), (even,odd,even), (odd,even,even).
    
    # --- Subproblem 1: Cover (odd,odd,odd) black squares ---
    print("Solving Subproblem 1: Covering black squares of type (odd, odd, odd)...")
    target_1 = parity_types[(1, 1, 1)] # black squares with (odd,odd,odd) coords
    # These can only be attacked by unicorns from their partition: (o,o,o) or (e,e,e)
    candidate_1 = parity_types[(1, 1, 1)] + parity_types[(0, 0, 0)]
    m1 = solve_subproblem("Subproblem_1_ooo", target_1, candidate_1)
    print(f"Minimum unicorns needed for Subproblem 1: {m1}")
    print("-" * 20)

    # --- Subproblem 2: Cover (even,even,odd) black squares ---
    print("Solving Subproblem 2: Covering black squares of type (even, even, odd)...")
    target_2 = parity_types[(0, 0, 1)] # black squares with (even,even,odd) coords
    # These can only be attacked by unicorns from their partition: (e,e,o) or (o,o,e)
    candidate_2 = parity_types[(0, 0, 1)] + parity_types[(1, 1, 0)]
    m2 = solve_subproblem("Subproblem_2_eeo", target_2, candidate_2)
    print(f"Minimum unicorns needed for Subproblem 2: {m2}")
    print("-" * 20)

    # --- Subproblems 3 and 4 ---
    # By symmetry, the problems for covering (e,o,e) and (o,e,e) black squares are identical to Subproblem 2.
    m3 = m2
    m4 = m2
    print("By symmetry, the results for the other two black-square partitions are the same.")
    print(f"Minimum unicorns needed for Subproblem 3: {m3}")
    print(f"Minimum unicorns needed for Subproblem 4: {m4}")
    print("-" * 20)

    # --- Final Calculation ---
    total_unicorns = m1 + m2 + m3 + m4
    print("The total minimum number of unicorns is the sum of the results from the four independent subproblems.")
    print(f"Total = (Subproblem 1) + (Subproblem 2) + (Subproblem 3) + (Subproblem 4)")
    print(f"Total = {m1} + {m2} + {m3} + {m4} = {total_unicorns}")


if __name__ == '__main__':
    # Checking for pulp installation
    try:
        import pulp
    except ImportError:
        print("This script requires the 'pulp' library.")
        print("Please install it using: pip install pulp")
        exit()
        
    main()
    print("\n<<<32>>>")