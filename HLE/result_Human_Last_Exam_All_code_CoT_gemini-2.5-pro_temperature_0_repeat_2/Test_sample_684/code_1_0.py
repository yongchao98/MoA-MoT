def find_solvable_hypercube_sizes():
    """
    This function finds and prints the integer values of n in the range [7, 55]
    for which it's possible to leave just one gift on the hypercube.

    The condition for this to be possible is that n must be congruent to 1 modulo 7.
    """
    
    # The range for n is 7 <= n <= 55.
    min_n = 7
    max_n = 55
    
    solvable_n_values = []
    
    # Iterate through all possible values of n.
    for n in range(min_n, max_n + 1):
        # The problem is solvable if and only if n % 7 == 1.
        if n % 7 == 1:
            solvable_n_values.append(n)
            
    # Print the values in increasing order.
    # The values are naturally found in increasing order.
    print("The values of n for which it is possible to reach a state with only one gift are:")
    print(*solvable_n_values, sep="\n")

find_solvable_hypercube_sizes()