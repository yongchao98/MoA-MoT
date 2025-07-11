def find_possible_n_values():
    """
    This function finds and prints the values of n for which it's possible
    to leave exactly one gift on the 5D hypercube.
    
    The problem is solvable if and only if n is of the form 14k + 1.
    We need to find such values of n in the range [7, 55].
    """
    
    # List to store the valid n values
    valid_n_values = []
    
    # Iterate through the given range for n
    for n in range(7, 56):
        # A solution exists if and only if n can be expressed as n = 14k + 1 for some integer k.
        # This is equivalent to n % 14 == 1.
        # This condition is derived from a linear algebra analysis over F_2, where
        # we require q = n // 7 to be even and r = n % 7 to be 1.
        if n % 14 == 1:
            valid_n_values.append(n)
            
    # Print each valid value of n. The loop guarantees they are in increasing order.
    for val in valid_n_values:
        print(val)

# Execute the function to find and print the answers.
find_possible_n_values()