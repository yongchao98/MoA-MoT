import collections

def solve_permutation_and_modulo():
    """
    This function finds the optimal permutation of a list of numbers 'a'
    to maximize the final value of 'x' after a series of modulo operations.
    It uses dynamic programming on subsets (masks) to explore all possible outcomes.
    """
    # Example Input. You can change these values.
    a = [6, 7, 8]
    x = 100
    
    n = len(a)
    
    # dp[mask] will be a dictionary {value: (previous_value, last_index)}
    # This allows reconstructing the path.
    dp = collections.defaultdict(dict)
    dp[0] = {x: (None, None)}

    for mask in range(1, 1 << n):
        for i in range(n):
            # If i-th element is in the current subset (mask)
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                # For each possible value from the previous state
                for val, _ in dp[prev_mask].items():
                    new_val = val % a[i]
                    # We store the path to this new_val. If multiple paths
                    # lead to the same new_val, we just keep one.
                    dp[mask][new_val] = (val, i)

    full_mask = (1 << n) - 1
    if not dp[full_mask]:
        # This case happens if a contains 0 or other issues,
        # though problem constraints usually prevent this.
        if n == 0:
            print(f"{x} = {x}")
            return
        else:
             print("No valid permutation found.")
             return


    # Find the best final value (the maximum one)
    best_val = -1
    for val in dp[full_mask].keys():
        if val > best_val:
            best_val = val

    # Reconstruct the permutation
    p = []
    current_val = best_val
    current_mask = full_mask
    for _ in range(n):
        prev_val, last_idx = dp[current_mask][current_val]
        p.append(a[last_idx])
        current_val = prev_val
        current_mask ^= (1 << last_idx)
    
    p.reverse()

    # Print the resulting equation
    equation_str = f"({x}"
    temp_x = x
    for val in p:
        temp_x %= val
        equation_str += f" % {val}"
    equation_str += f") = {temp_x}"

    print(f"To minimize the absolute difference |{x} - final_x|, we need to maximize final_x.")
    print(f"The best permutation found is: {p}")
    print("The equation is:")
    
    # Print the full equation with intermediate steps for clarity
    print_equation = f"{x}"
    final_result_check = x
    for num in p:
        print_equation += f" % {num}"
        final_result_check %= num
    print_equation += f" = {final_result_check}"
    
    print(print_equation)


solve_permutation_and_modulo()
