import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Generate all valid row masks for a row of length m.
    # A mask is valid if no two adjacent bits are 1.
    valid_masks = []
    for i in range(1 << m):
        if (i & (i << 1)) == 0:
            valid_masks.append(i)
    
    num_valid_masks = len(valid_masks)

    # Step 2: Create the transition matrix T.
    # T[i][j] = 1 if mask i and mask j can be in adjacent rows.
    # Let T[i][j] mean previous row is mask_i, current row is mask_j
    T = np.zeros((num_valid_masks, num_valid_masks), dtype=np.uint64)
    mask_to_idx = {mask: i for i, mask in enumerate(valid_masks)}

    for i, prev_mask in enumerate(valid_masks):
        for j, current_mask in enumerate(valid_masks):
            # The vertical constraint: no 1 above a 1.
            if (prev_mask & current_mask) == 0:
                T[i, j] = 1

    # Step 3: Dynamic Programming Calculation
    # counts[k] stores the number of ways to fill rows up to the current one,
    # ending with valid_masks[k].
    
    # Step 4: Base case for the first row.
    # Any valid mask can be the first row.
    counts = np.ones(num_valid_masks, dtype=np.uint64)
    
    # Iterate for the remaining n-1 rows
    for _ in range(n - 1):
        # The number of ways for the current row with mask j is the sum of ways
        # for the previous row over all compatible masks i.
        # This is a matrix-vector multiplication: counts = counts * T
        counts = counts.dot(T)

    # Step 5: The final 'counts' array contains the number of ways to form a valid
    # n x m grid, for each possible valid last row. The total is their sum.
    total_assignments = np.sum(counts)

    final_sum_str = " + ".join(map(str, counts))
    print(f"The number of valid assignments for a {n}x{m} grid is found by summing the number of ways to form a valid grid ending in each possible configuration for the last row.")
    print("\nCalculation:")
    print(f"{final_sum_str} = {total_assignments}")

solve_grid_assignments()