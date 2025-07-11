import math

def solve():
    """
    Finds the permutation of a list of numbers 'a' that minimizes the absolute 
    difference between an initial value 'x' and the final value after applying 
    a sequence of modulo operations.
    """
    
    # Example Input
    a = [10, 5, 17]
    x_initial = 25
    
    n = len(a)

    if n > 20:
        print("Input list is too large (n > 20). This algorithm is exponential and may be too slow.")
        return

    # dp[mask] = {value: (prev_value, index_of_number_used)}
    # This structure allows us to backtrack the path.
    dp = {}
    dp[0] = {x_initial: (None, None)}

    for mask in range(1, 1 << n):
        dp[mask] = {}
        for i in range(n):
            # Check if the i-th number is in the current subset (mask)
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                # For each value in the previous state, calculate the new value
                if prev_mask in dp:
                    for prev_val, _ in dp[prev_mask].items():
                        new_val = prev_val % a[i]
                        # Store the new value and how we got here
                        dp[mask][new_val] = (prev_val, i)

    final_mask = (1 << n) - 1
    if not dp.get(final_mask):
        print("No solution found. This can happen if the input list is empty.")
        return
        
    final_values = dp[final_mask].keys()

    # Find the final value that is closest to the initial x
    best_x = min(final_values, key=lambda val: abs(val - x_initial))

    # Backtrack to find the permutation
    p = []
    current_x = best_x
    current_mask = final_mask

    for _ in range(n):
        prev_val, index = dp[current_mask][current_x]
        p.append(a[index])
        current_x = prev_val
        current_mask ^= (1 << index)
    
    # The permutation is built backwards, so we reverse it
    p.reverse()

    # Print the equation
    equation_str = f"{x_initial} % " + " % ".join(map(str, p))
    print(f"{equation_str} = {best_x}")

solve()