import collections

def solve():
    """
    This function solves the permutation modulo problem using dynamic programming.
    It finds the permutation of 'a' that maximizes the final value of 'x'
    after sequential modulo operations.
    """
    # Example Input
    a = [3, 5, 8]
    x = 20

    n = len(a)
    full_mask = (1 << n) - 1

    # dp[mask] will store the set of reachable values using numbers corresponding to the mask
    dp = collections.defaultdict(set)
    dp[0] = {x}

    # Fill the DP table
    for mask in range(1, 1 << n):
        for i in range(n):
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                if prev_mask in dp:
                    for val in dp[prev_mask]:
                        dp[mask].add(val % a[i])

    # Find the best final x
    if not dp[full_mask]:
        best_x = x # Case where list 'a' is empty
    else:
        best_x = max(dp[full_mask])

    # Reconstruct one of the permutations that yields best_x
    permutation = []
    current_x = best_x
    current_mask = full_mask

    for _ in range(n, 0, -1):
        found_prev = False
        for i in range(n):
            if (current_mask >> i) & 1:
                prev_mask = current_mask ^ (1 << i)
                # Check if some value in the previous state can lead to current_x
                if prev_mask in dp:
                    for prev_x in dp[prev_mask]:
                        if prev_x % a[i] == current_x:
                            permutation.append(a[i])
                            current_x = prev_x
                            current_mask = prev_mask
                            found_prev = True
                            break
            if found_prev:
                break
    
    permutation.reverse()

    # Print the equation step by step
    temp_x = x
    equation_parts = [str(x)]
    for val in permutation:
        next_x = temp_x % val
        equation_parts.append(f"mod {val} -> {next_x}")
        temp_x = next_x
        
    final_equation = " ".join(equation_parts)
    print(f"Initial list: a = {a}")
    print(f"Initial x = {x}")
    print(f"An optimal permutation found: {permutation}")
    print(f"Calculation: {final_equation}")
    print(f"The best resulting x is: {best_x}")
    print(f"The smallest absolute difference is: {abs(x - best_x)}")


solve()