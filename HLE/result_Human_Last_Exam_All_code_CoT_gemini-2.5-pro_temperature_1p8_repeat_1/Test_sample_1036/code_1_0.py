import math

# Use a cache for gcd to speed up calculations
gcd_cache = {}
def gcd(a, b):
    if (a, b) in gcd_cache:
        return gcd_cache[(a, b)]
    res = math.gcd(a, b)
    gcd_cache[(a, b)] = res
    return res

# The recursive backtracking function
def count_valid_sequences(k, remainders):
    """
    Counts the number of valid sequences of remainders recursively.
    k: the current divisor to choose a remainder for.
    remainders: a dictionary {divisor: remainder} for previously chosen values.
    """
    if k > 100:
        # Base case: successfully found a valid sequence up to 100
        return 1

    count = 0
    # Iterate through all possible values for the remainder r_k
    for r_k in range(k):
        # Check for distinctness from previously chosen remainders
        if r_k in remainders.values():
            continue

        is_consistent = True
        # Check for consistency with all previous remainders
        for j, r_j in remainders.items():
            if r_k % gcd(k, j) != r_j:
                is_consistent = False
                break
        
        if is_consistent:
            # If r_k is valid, proceed to the next divisor k+1
            new_remainders = remainders.copy()
            new_remainders[k] = r_k
            count += count_valid_sequences(k + 1, new_remainders)
            
    return count

# Start the calculation
result = count_valid_sequences(2, {})

# Print the final result
print("The number of such positive integers is:")
print(result)
