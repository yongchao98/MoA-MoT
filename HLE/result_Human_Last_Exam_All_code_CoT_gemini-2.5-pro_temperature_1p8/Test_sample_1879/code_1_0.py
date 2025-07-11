import math

def combinations(n, k):
    """Helper function to calculate combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_grid_path():
    """
    Calculates the number of unique ways to move from (0,0) to (4,8)
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_steps = 4
    U_steps = 8
    total_steps = R_steps + U_steps
    
    print("This problem can be solved using combinatorics and the Principle of Inclusion-Exclusion.")
    print(f"We are moving from (0,0) to ({R_steps},{U_steps}), which requires {R_steps} Right moves and {U_steps} Up moves.\n")

    # Step 1: Calculate total paths without any restrictions
    total_paths = combinations(total_steps, R_steps)
    print("Step 1: Total number of paths without restrictions.")
    print(f"This is the number of ways to arrange {R_steps} 'R's in {total_steps} total steps.")
    print(f"Total Paths = C({total_steps}, {R_steps}) = {total_paths}\n")

    # Step 2: Calculate the number of invalid paths (containing 'RRRR' or 'UUUU')
    print("Step 2: Calculate the number of invalid paths.")
    print("An invalid path contains 'RRRR' (Set A) or 'UUUU' (Set B).")
    print("We need to find |A U B| = |A| + |B| - |A intersect B|.\n")
    
    # Step 2a: Calculate |A| (paths with 'RRRR')
    # Treat 'RRRR' as a single block. We arrange this block with 8 'U's.
    # Total items to arrange = 1 (RRRR block) + 8 ('U's) = 9
    # Choose 1 position for the RRRR block out of 9.
    paths_with_rrrr = combinations(U_steps + 1, 1)
    print("Step 2a: Calculate |A| (paths with 'RRRR').")
    print("Treat 'RRRR' as a single block. We arrange this block and the 8 'U's.")
    print(f"|A| = C({U_steps} + 1, 1) = C(9, 1) = {paths_with_rrrr}\n")

    # Step 2b: Calculate |B| (paths with 'UUUU')
    # Use stars and bars. The 4 R's are bars, creating 5 slots for 8 U's (stars).
    # u1+u2+u3+u4+u5 = 8. A path is invalid if any ui >= 4.
    # Use inclusion-exclusion on the slots.
    slots = R_steps + 1
    # Term 1: At least one slot has >= 4 'U's
    term1_comb_n = U_steps - 4 + slots - 1
    term1 = combinations(slots, 1) * combinations(term1_comb_n, slots - 1)
    # Term 2: At least two slots have >= 4 'U's (i.e., 4+4=8)
    term2_comb_n = U_steps - 8 + slots - 1
    term2 = combinations(slots, 2) * combinations(term2_comb_n, slots - 1)
    
    paths_with_uuuu = term1 - term2
    
    print("Step 2b: Calculate |B| (paths with 'UUUU').")
    print("We use stars and bars, where the 4 'R's are bars and 8 'U's are stars.")
    print("This creates 5 slots for the 'U's. We count configurations where any slot has >= 4 'U's.")
    print(f"Using Inclusion-Exclusion for slots:")
    print(f"  Ways with >=1 slot violating = C(5,1) * C(8-4 + 5-1, 5-1) = {combinations(slots, 1)} * {combinations(term1_comb_n, slots - 1)} = {term1}")
    print(f"  Ways with >=2 slots violating = C(5,2) * C(8-8 + 5-1, 5-1) = {combinations(slots, 2)} * {combinations(term2_comb_n, slots - 1)} = {term2}")
    print(f"|B| = {term1} - {term2} = {paths_with_uuuu}\n")
    
    # Step 2c: Analyze the intersection
    # A path with RRRR places 8 U's in 2 slots (before and after RRRR). u1+u2=8.
    # For this path *not* to have UUUU, both u1<=3 and u2<=3, which is impossible (3+3=6!=8).
    # So every path in A is also in B. A is a subset of B.
    # Therefore |A U B| = |B|.
    invalid_paths = paths_with_uuuu
    print("Step 2c: Calculate the total number of invalid paths.")
    print("Any path with 'RRRR' must place the 8 'U's around it. This forces a run of at least 5 'U's.")
    print("This means any path in Set A is also in Set B (A is a subset of B).")
    print(f"Therefore, the total number of invalid paths is |A U B| = |B| = {invalid_paths}\n")
    
    # Step 3: Final calculation
    valid_paths = total_paths - invalid_paths
    print("Step 3: Calculate the final answer.")
    print("Valid Paths = Total Paths - Invalid Paths")
    print(f"Final Answer = {total_paths} - {invalid_paths} = {valid_paths}")
    
solve_grid_path()
<<<155>>>