import sys
from collections import Counter

def solve():
    """
    Solves the modulo permutation problem.
    """
    # Example Input:
    # Given a list of n numbers a_1, a_2, ..., a_n and another number x
    a = [3, 5, 8]
    x_orig = 10

    print(f"The given list is: {a}")
    print(f"The initial value of x is: {x_orig}")

    if not a:
        print("\nThe list is empty.")
        print(f"Final x equals initial x: {x_orig}")
        print(f"Equation: {x_orig}")
        print(f"Absolute difference: 0")
        return

    # 1. Simplify the problem by finding the minimum element 'm'
    m = min(a)

    # Create a_prime by removing one instance of m.
    # We use a Counter to handle duplicates correctly.
    a_counts = Counter(a)
    a_prime = []
    a_counts[m] -= 1
    for num, count in a_counts.items():
        a_prime.extend([num] * count)
    
    n_prime = len(a_prime)

    # 2. Dynamic Programming on subsets of a_prime
    # dp[mask] will store the set of reachable values using the subset 'mask'
    dp = [set() for _ in range(1 << n_prime)]
    dp[0] = {x_orig}

    for mask in range(1, 1 << n_prime):
        for i in range(n_prime):
            # Check if the i-th element of a_prime is in the current subset
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                mod_val = a_prime[i]
                for val in dp[prev_mask]:
                    dp[mask].add(val % mod_val)

    # 3. Collect all possible intermediate values before applying 'm'
    all_intermediate_values = set()
    for mask_set in dp:
        all_intermediate_values.update(mask_set)

    # Calculate final values by taking modulo m. If a_prime is empty, the only intermediate is x_orig.
    if not all_intermediate_values:
        final_values = {x_orig % m}
    else:
        final_values = {val % m for val in all_intermediate_values}

    # 4. Find the best final x that minimizes the absolute difference
    best_x = -1
    min_diff = float('inf')

    for val in sorted(list(final_values)) :
        diff = abs(x_orig - val)
        if diff < min_diff:
            min_diff = diff
            best_x = val

    # 5. Reconstruct the permutation for the equation
    # Find an intermediate value `v` and a `mask` that produced `best_x`
    source_v = -1
    source_mask = -1
    
    # We search for a `v` that results in `best_x`.
    # To make reconstruction deterministic, we can sort the values.
    for v in sorted(list(all_intermediate_values), reverse=True):
        if v % m == best_x:
            # Find the smallest mask that can generate this `v`
            for i in range(1 << n_prime):
                if v in dp[i]:
                    source_v = v
                    source_mask = i
                    break
            if source_v != -1:
                break

    # If no intermediate values (e.g., a_prime is empty), handle separately
    if source_v == -1:
         path = [m]
         source_mask=0
    else:
        path = [m]
        curr_v = source_v
        curr_mask = source_mask

        # Backtrack through the DP table
        while curr_mask > 0:
            found_prev = False
            for i in range(n_prime):
                if (curr_mask >> i) & 1:
                    prev_mask = curr_mask ^ (1 << i)
                    mod_val = a_prime[i]
                    # We need to find the predecessor value
                    for prev_v in dp[prev_mask]:
                        if prev_v % mod_val == curr_v:
                            path.append(mod_val)
                            curr_v = prev_v
                            curr_mask = prev_mask
                            found_prev = True
                            break
                    if found_prev:
                        break
            if not found_prev:
                break # Should not happen in a correct trace
        path.reverse()

    # Construct the full permutation
    path_counts = Counter(path)
    full_permutation = list(path)
    # Add unused numbers to the permutation (their order doesn't matter)
    for num, count in Counter(a).items():
        used_count = path_counts.get(num, 0)
        if count > used_count:
            full_permutation.extend([num] * (count - used_count))

    # Print the equation
    equation = f"{x_orig}"
    temp_x = x_orig
    for num in full_permutation:
        equation += f" mod {num}"
        temp_x %= num
    equation += f" = {temp_x}"

    print(f"\nThe best resulting x is: {best_x}")
    print(f"The minimum absolute difference is: {min_diff}")
    print(f"An optimal permutation is: {full_permutation}")
    print(f"The corresponding equation is: {equation}")

solve()