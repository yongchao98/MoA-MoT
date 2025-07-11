def solve_modulo_permutation():
    """
    Finds the permutation of a list 'a' that maximizes the final value of 'x'
    after sequential a series of modulo operations. It uses dynamic programming
    on subsets (bitmask DP).
    """
    # Example input. You can change these values to test with other cases.
    x_initial = 100
    a = [11, 21, 31]
    n = len(a)

    if n == 0:
        print(f"{x_initial} = {x_initial}")
        return
        
    # dp[mask] stores the maximum final x for the subset of numbers represented by mask.
    dp = {}
    # path[mask] stores the index of the last element in the optimal sequence for mask.
    path = {}

    # Base cases: masks with a single element (subset size 1)
    for i in range(n):
        mask = 1 << i
        dp[mask] = x_initial % a[i]
        path[mask] = i

    # Iterate through subset sizes from 2 to n
    for size in range(2, n + 1):
        # Iterate through all subsets of the current size
        for mask in range(1, 1 << n):
            if bin(mask).count('1') == size:
                dp[mask] = -1
                # For each subset, try every element 'a[j]' as the last one in the sequence
                for j in range(n):
                    if (mask >> j) & 1:
                        prev_mask = mask ^ (1 << j)
                        
                        # The value before the last operation is the optimal result
                        # for the subproblem on the set represented by prev_mask.
                        val = dp[prev_mask] % a[j]
                        
                        if val > dp[mask]:
                            dp[mask] = val
                            path[mask] = j
    
    final_mask = (1 << n) - 1
    x_final = dp[final_mask]

    # Reconstruct the permutation by backtracking from the final state
    p_indices = []
    curr_mask = final_mask
    while curr_mask != 0:
        last_element_index = path[curr_mask]
        p_indices.append(last_element_index)
        curr_mask ^= (1 << last_element_index)
    
    # The permutation is found in reverse order, so we reverse it back.
    p_indices.reverse()
    
    # Construct and print the final equation string showing each step
    equation_parts = [str(x_initial)]
    for index in p_indices:
        equation_parts.append(f"% {a[index]}")
    
    equation = " ".join(equation_parts)
    equation += f" = {x_final}"
    
    print(equation)

solve_modulo_permutation()