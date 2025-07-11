def solve():
    """
    Finds all values of n between 7 and 55 (inclusive) for which it's possible
    to have only one gift remaining on an n x n x n x n x n hypercube.

    The condition for this to be possible is that n modulo 7 must be 6.
    """
    
    # The list to store the valid n values
    valid_n_values = []
    
    # Iterate through the given range for n
    for n in range(7, 56):
        # Check if n is congruent to 6 modulo 7
        if n % 7 == 6:
            valid_n_values.append(n)
            
    # Print the result
    print("The possible values for n are:")
    print(valid_n_values)

solve()