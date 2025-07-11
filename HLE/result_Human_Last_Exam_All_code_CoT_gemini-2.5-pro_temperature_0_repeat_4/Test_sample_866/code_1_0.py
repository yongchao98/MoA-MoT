import math

def solve_flat_folding():
    """
    Calculates the total number of different assignments of creases 
    which will make the given pattern flat fordable.
    """
    # The input pattern of angles and creases
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
    
    # Step 1: Extract the crease assignments from the pattern
    creases = [p for p in pattern if isinstance(p, str)]
    print(f"The sequence of creases is: {creases}")
    
    # Step 2: Count the number of known Mountain, Valley, and Unassigned folds
    M_known = creases.count('M')
    V_known = creases.count('V')
    n_q = creases.count('?')
    
    print(f"Initial counts: {M_known} Mountain, {V_known} Valley, {n_q} Unassigned.")
    print("-" * 50)
    
    total_valid_assignments = 0
    
    # Step 3: Iterate through all possible ways to assign the '?' folds.
    # 'k' represents the number of '?'s that are assigned as Mountain folds.
    for k in range(n_q + 1):
        num_new_M = k
        num_new_V = n_q - k
        
        # Step 4: Calculate the total number of M and V folds for this assignment
        M_total = M_known + num_new_M
        V_total = V_known + num_new_V
        
        # The Maekawa-Kawasaki theorem for a flat-foldable vertex is |M - V| = 2
        diff = abs(M_total - V_total)
        
        print(f"Checking assignment of {num_new_M} '?' as 'M' and {num_new_V} '?' as 'V':")
        # Output the numbers in the final equation as requested
        print(f"  Equation check: |M_total - V_total| = |{M_total} - {V_total}| = {diff}")
        
        # Step 5: Check if the flat-foldability condition is met
        if diff == 2:
            # This assignment is valid. Calculate the number of ways it can be done.
            # This is C(n_q, k), i.e., "n_q choose k".
            num_ways = math.comb(n_q, k)
            print(f"  Result: Condition MET. Number of ways for this case: C({n_q}, {k}) = {num_ways}")
            total_valid_assignments += num_ways
        else:
            print(f"  Result: Condition NOT met.")
        print("-" * 50)

    # Step 6: Print the final total count
    print(f"The total number of different flat-foldable assignments is: {total_valid_assignments}")

solve_flat_folding()
<<<0>>>