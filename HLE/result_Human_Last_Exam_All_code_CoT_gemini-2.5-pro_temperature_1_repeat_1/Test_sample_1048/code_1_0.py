import sys

def solve():
    """
    Solves the modulo permutation problem using dynamic programming on subsets.
    """
    # Read input from the user.
    # Example format: x = 100, a = 10 20 30
    try:
        line = input("Enter x and the list of numbers a (e.g., x = 100, a = 10 20 30): ")
        x_str, a_str = line.split(',')
        initial_x = int(x_str.split('=')[1].strip())
        a = [int(num) for num in a_str.split('=')[1].strip().split()]
        n = len(a)
        
        if n == 0:
            print(f"The list 'a' is empty. The final x is {initial_x}.")
            return
            
    except (ValueError, IndexError):
        print("Invalid input format. Please use the format 'x = <value>, a = <values separated by space>'.")
        print("Using default example: x = 50, a = [15, 25, 35]")
        initial_x = 50
        a = [15, 25, 35]
        n = len(a)

    # dp[mask] will be a dictionary {value: (previous_value, last_element_index)}
    # This stores how each value was reached.
    dp = [{} for _ in range(1 << n)]
    dp[0] = {initial_x: (None, None)}

    # Fill the DP table
    for mask in range(1, 1 << n):
        for i in range(n):
            # Check if the i-th element is in the current subset
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                if not dp[prev_mask]:
                    continue
                
                # For each value reachable with the previous subset,
                # apply the modulo operation with a[i].
                for val, prev_state in dp[prev_mask].items():
                    new_val = val % a[i]
                    # Store one way to reach new_val for this mask
                    dp[mask][new_val] = (val, i)

    full_mask = (1 << n) - 1
    if not dp[full_mask]:
        print("No solution found. This can happen if the list 'a' is empty.")
        # This case is handled above, but as a safeguard.
        return
        
    # Find the best final value
    best_final_x = max(dp[full_mask].keys())

    # Reconstruct the permutation
    p = []
    curr_val = best_final_x
    curr_mask = full_mask
    while curr_mask > 0:
        prev_val, prev_idx = dp[curr_mask][curr_val]
        p.append(a[prev_idx])
        curr_val = prev_val
        curr_mask ^= (1 << prev_idx)
    
    p.reverse()

    # Build the equation string
    equation_str = f"{best_final_x} = "
    temp_x_str = str(initial_x)
    for val in p:
        temp_x_str = f"({temp_x_str} % {val})"
    equation_str += temp_x_str

    print(f"Original x: {initial_x}")
    print(f"List a: {a}")
    print(f"Optimal permutation: {p}")
    print(f"Final x: {best_final_x}")
    print("Equation:")
    print(equation_str)

solve()