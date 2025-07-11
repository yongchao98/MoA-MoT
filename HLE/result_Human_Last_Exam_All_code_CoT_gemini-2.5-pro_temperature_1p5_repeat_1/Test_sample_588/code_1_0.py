def get_chi(a, b):
    """Calculates the Euler characteristic of M(a) x M(b)."""
    return (2 - 2 * a) * (2 - 2 * b)

def find_solution():
    """
    Finds the lexicographically least tuple (a1, b1, ..., al, bl) 
    with l minimal, such that no M(ai, bi) is full, yet their 
    connect-sum is full.
    """
    # We determined that the minimal length l is 2.
    # We need to find two pairs (a1, b1) and (a2, b2)
    # such that a,b are not 1, and chi(a1,b1) + chi(a2,b2) = 0.
    
    search_limit = 5 # A small limit is sufficient for the smallest pairs.

    # Generate candidate pairs (a,b) with a,b != 1 and a <= b, sorted lexicographically.
    candidate_pairs = []
    for a in range(search_limit):
        if a == 1:
            continue
        for b in range(a, search_limit):
            if b == 1:
                continue
            candidate_pairs.append((a, b))

    # Iterate through pairs to find the first (lexicographically smallest) solution
    for p1 in candidate_pairs:
        chi1 = get_chi(p1[0], p1[1])
        target_chi = -chi1
        
        for p2 in candidate_pairs:
            if get_chi(p2[0], p2[1]) == target_chi:
                # Solution found. Order the pairs and form the final tuple.
                solution_pairs = sorted([p1, p2])
                
                final_tuple = solution_pairs[0] + solution_pairs[1]
                
                print("Found two manifolds M(a,b) which are not full, but whose connect-sum is full.")
                print(f"1. M({solution_pairs[0][0]},{solution_pairs[0][1]}): Euler characteristic = {get_chi(solution_pairs[0][0], solution_pairs[0][1])}")
                print(f"2. M({solution_pairs[1][0]},{solution_pairs[1][1]}): Euler characteristic = {get_chi(solution_pairs[1][0], solution_pairs[1][1])}")
                
                chi_val1 = get_chi(solution_pairs[0][0], solution_pairs[0][1])
                chi_val2 = get_chi(solution_pairs[1][0], solution_pairs[1][1])
                
                print(f"\nThe final equation for the Euler characteristic of the connect-sum is:")
                print(f"{chi_val1} + ({chi_val2}) = {chi_val1 + chi_val2}")
                
                print("\nThe lexicographically least tuple is:")
                # Output the tuple in the required flat format
                print(f"({','.join(map(str, final_tuple))})")
                return

find_solution()
