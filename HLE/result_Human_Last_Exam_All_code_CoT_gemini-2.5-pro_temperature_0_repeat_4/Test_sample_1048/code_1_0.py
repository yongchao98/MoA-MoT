import sys

def solve():
    """
    Solves the permutation of remainders problem using dynamic programming on subsets.
    """
    try:
        # Read inputs from the user
        n_str = input("Enter the number of elements in the list (n): ")
        n = int(n_str)
        
        a_str = input("Enter the list of n numbers, separated by spaces: ")
        a = [int(num) for num in a_str.split()]
        
        x_str = input("Enter the initial number x: ")
        x = int(x_str)

        if len(a) != n:
            print(f"Error: Expected {n} numbers, but got {len(a)}.", file=sys.stderr)
            return

    except (ValueError, IndexError):
        print("Invalid input. Please enter valid integers.", file=sys.stderr)
        return

    # Sort a to optimize the DP state size and keep original indices for the permutation
    # The problem asks for a permutation of the list, so original indices are not strictly
    # necessary, but sorting is key for the algorithm's performance.
    a.sort()

    # dp[mask] = set of achievable remainders using numbers in the mask
    dp = [set() for _ in range(1 << n)]
    dp[0] = {x}

    # parent[mask][new_val] = (index_of_a_i, prev_val) to reconstruct the path
    parent = [{} for _ in range(1 << n)]

    for mask in range(1, 1 << n):
        for i in range(n):
            if (mask >> i) & 1:  # If i-th element is in the current subset
                prev_mask = mask ^ (1 << i)
                if not dp[prev_mask]:
                    continue
                
                a_i = a[i]
                for prev_val in dp[prev_mask]:
                    new_val = prev_val % a_i
                    dp[mask].add(new_val)
                    # Store one possible path. Since we need just one permutation, this is fine.
                    parent[mask][new_val] = (i, prev_val)

    final_mask = (1 << n) - 1
    if not dp[final_mask]:
        print("No solution found. This can happen if the initial list is empty.")
        # If n=0, the loop doesn't run. The final x is the original x.
        if n == 0:
            print(f"Initial x = {x}")
            print("Optimal permutation of a: []")
            print(f"Final x = {x}")
        return

    # Find the best resulting x (the maximum one)
    best_x = -1
    # Using max() on an empty set raises an error, so handle it.
    if dp[final_mask]:
        best_x = max(dp[final_mask])

    # Reconstruct the permutation by backtracking
    p = []
    curr_mask = final_mask
    curr_val = best_x
    
    # We need to backtrack n steps
    for _ in range(n):
        if curr_mask == 0: break
        i, prev_val = parent[curr_mask][curr_val]
        p.append(a[i])
        curr_mask ^= (1 << i)
        curr_val = prev_val
    
    p.reverse()

    # Print the results as requested
    print("\n--- Calculation ---")
    print(f"Initial x = {x}")
    print(f"Optimal permutation of a: {p}")
    
    # Build the equation string and print step-by-step calculation
    final_eq_str = f"{x}"
    current_x_val = x
    
    print("\nStep-by-step operations:")
    for val in p:
        next_x_val = current_x_val % val
        print(f"{current_x_val} mod {val} = {next_x_val}")
        final_eq_str = f"({final_eq_str}) mod {val}"
        current_x_val = next_x_val
        
    print("\n--- Final Result ---")
    print("Final equation:")
    print(f"{final_eq_str} = {best_x}")
    print(f"The resulting x is: {best_x}")
    print(f"The absolute difference |{best_x} - {x}| is {abs(best_x - x)}")


solve()