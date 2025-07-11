import collections

memo_comp = {}

def count_compositions(n, k, max_val):
    """
    Calculates the number of ways to write n as a sum of k ordered integers,
    each between 1 and max_val. This is equivalent to finding the number
    of compositions of n into k parts, with each part at most max_val.
    Uses memoization to avoid recomputing results.
    """
    if (n, k, max_val) in memo_comp:
        return memo_comp[(n, k, max_val)]
    
    # Base cases
    if k == 0:
        return 1 if n == 0 else 0
    # If the remaining sum 'n' is impossible to achieve with 'k' parts
    if n < k or k * max_val < n:
        return 0
    
    # Recursive step: sum the possibilities for the last part's size
    count = 0
    for part_size in range(1, max_val + 1):
        if n - part_size >= 0:
            count += count_compositions(n - part_size, k - 1, max_val)
            
    memo_comp[(n, k, max_val)] = count
    return count

def solve_path_problem():
    """
    Calculates the number of valid paths from (0,0) to (4,8).
    """
    total_rights = 4
    total_ups = 8
    max_run = 3
    
    max_k = total_rights + total_ups

    R = collections.defaultdict(int)
    U = collections.defaultdict(int)

    for k in range(1, max_k + 1):
        if k <= total_rights:
            R[k] = count_compositions(total_rights, k, max_run)
        if k <= total_ups:
            U[k] = count_compositions(total_ups, k, max_run)

    # Case 1: k R-blocks and k U-blocks
    term1 = 0
    for k in range(1, max_k + 1):
        term1 += 2 * R[k] * U[k]

    # Case 2: k R-blocks and k-1 U-blocks
    term2 = 0
    for k in range(1, max_k + 1):
        term2 += R[k] * U[k-1]

    # Case 3: k U-blocks and k-1 R-blocks
    term3 = 0
    for k in range(1, max_k + 1):
        term3 += U[k] * R[k-1]

    total_ways = term1 + term2 + term3
    
    print(f"The total number of unique ways is the sum of three cases:")
    print(f"1. Paths with equal numbers of R-blocks and U-blocks: {term1}")
    print(f"2. Paths with more R-blocks than U-blocks: {term2}")
    print(f"3. Paths with more U-blocks than R-blocks: {term3}")
    print("\nFinal calculation:")
    print(f"{term1} + {term2} + {term3} = {total_ways}")
    print(f"\nTotal unique ways: {total_ways}")


solve_path_problem()