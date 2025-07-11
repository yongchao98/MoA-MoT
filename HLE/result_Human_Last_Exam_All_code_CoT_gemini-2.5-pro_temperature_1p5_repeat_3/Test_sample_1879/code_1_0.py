import math

def solve_grid_path_problem():
    """
    Calculates the number of unique paths from (0,0) to (4,8) with movement constraints.

    The method is to count the ways to partition the total 'Right' and 'Up' moves
    into runs of allowed lengths (1, 2, or 3) and then combine them based on how
    these runs can be arranged to form a valid path.
    """

    def combinations(n, k):
        """Calculates C(n, k), the number of combinations."""
        if k < 0 or k > n:
            return 0
        try:
            return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))
        except ValueError:
            return 0

    def count_compositions(n, k, max_val):
        """
        Counts compositions of n into k parts, where each part is >= 1 and <= max_val.
        This uses the principle of inclusion-exclusion.
        """
        total = 0
        for j in range(k + 1):
            # C(n-1, k-1) gives compositions with no max_val limit.
            # We subtract cases where one or more parts are > max_val.
            # n' = n - j*max_val effectively enforces this.
            term = ((-1)**j) * combinations(k, j) * combinations(n - j * max_val - 1, k - 1)
            total += term
        return total

    total_R, total_U = 4, 8
    max_run_size = 3  # Cannot have 4 or more consecutive steps.

    # Step 1: Calculate C_R(k), the number of ways to partition 4 'R' moves
    # into k runs, each of size at most 3.
    cr_vals = {k: count_compositions(total_R, k, max_run_size) for k in range(1, total_R + 1)}

    # Step 2: Calculate C_U(k), the number of ways to partition 8 'U' moves
    # into k runs, each of size at most 3.
    cu_vals = {k: count_compositions(total_U, k, max_run_size) for k in range(1, total_U + 1)}
    
    print("This script calculates the number of unique paths from A(0,0) to B(4,8).")
    print("Movement is restricted to 1 unit right (R) or up (U), with no more than 3 consecutive steps in the same direction.\n")
    print(f"Number of ways to partition 4 'R' moves into k runs (C_R(k)): {cr_vals}")
    print(f"Number of ways to partition 8 'U' moves into k runs (C_U(k)): {cu_vals}\n")

    # Step 3: Sum the possibilities based on the arrangement of runs.

    # Case 1: Number of R-runs (k) == Number of U-runs (k).
    # Path can be R...U or U...R. The factor of 2 accounts for this.
    case1_total = 0
    case1_eq_parts = []
    # k must be valid for both R and U compositions. Min k is 3 since ceil(8/3)=3.
    for k in range(3, 5): 
        if cr_vals.get(k, 0) > 0 and cu_vals.get(k, 0) > 0:
            term = cr_vals[k] * cu_vals[k]
            case1_total += term
            case1_eq_parts.append(f"{cr_vals[k]} * {cu_vals[k]}")
    case1_final_total = 2 * case1_total
    
    print("Case 1: Paths with an equal number of R-runs and U-runs (e.g., R...U or U...R).")
    print(f"Ways = 2 * (C_R(3)*C_U(3) + C_R(4)*C_U(4))")
    print(f"     = 2 * ({' + '.join(case1_eq_parts)})")
    print(f"     = 2 * ({case1_total}) = {case1_final_total}\n")

    # Case 2: Number of R-runs (k_r) == Number of U-runs (k_u) + 1.
    # Path must start and end with 'R'.
    case2_total = 0
    case2_eq_parts = []
    # Iterate through possible k_r and find matching k_u = k_r - 1.
    for k_r in range(2, 5): 
        k_u = k_r - 1
        if cr_vals.get(k_r, 0) > 0 and cu_vals.get(k_u, 0) > 0:
             term = cr_vals[k_r] * cu_vals[k_u]
             case2_total += term
             case2_eq_parts.append(f"{cr_vals[k_r]} * {cu_vals[k_u]}")
    
    print("Case 2: Paths with one more R-run than U-runs (must start and end with R).")
    print(f"Ways = C_R(4)*C_U(3)")
    print(f"     = {' + '.join(case2_eq_parts if case2_eq_parts else ['0'])}")
    print(f"     = {case2_total}\n")

    # Case 3: Number of U-runs (k_u) == Number of R-runs (k_r) + 1.
    # Path must start and end with 'U'.
    case3_total = 0
    case3_eq_parts = []
    # Iterate through possible k_u and find matching k_r = k_u - 1.
    for k_u in range(3, 9): 
        k_r = k_u - 1
        if cu_vals.get(k_u, 0) > 0 and cr_vals.get(k_r, 0) > 0:
            term = cu_vals[k_u] * cr_vals[k_r]
            case3_total += term
            # Show C_R * C_U for consistency
            case3_eq_parts.append(f"{cr_vals[k_r]} * {cu_vals[k_u]}")
            
    print("Case 3: Paths with one more U-run than R-runs (must start and end with U).")
    print(f"Ways = C_R(2)*C_U(3) + C_R(3)*C_U(4) + C_R(4)*C_U(5)")
    print(f"     = {' + '.join(case3_eq_parts if case3_eq_parts else ['0'])}")
    print(f"     = {9 + 57 + 30} = {case3_total}\n")

    # Final Total
    total_ways = case1_final_total + case2_total + case3_total
    print("Total unique ways = (Ways from Case 1) + (Ways from Case 2) + (Ways from Case 3)")
    print(f"Final Answer = {case1_final_total} + {case2_total} + {case3_total} = {total_ways}")
    
    return total_ways

final_answer = solve_grid_path_problem()
print(f"<<<{final_answer}>>>")