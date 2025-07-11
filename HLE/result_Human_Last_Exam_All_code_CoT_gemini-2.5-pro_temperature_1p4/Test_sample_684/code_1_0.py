def find_solvable_hypercube_sizes():
    """
    Finds and prints the integer values of n in the range [7, 55]
    for which it's possible to leave a single gift on the hypercube.
    
    The condition for solvability is that n must be congruent to 1 modulo 7.
    """
    
    # The given range for n
    min_n = 7
    max_n = 55
    
    # List to store the valid values of n
    solvable_n_values = []
    
    # Iterate through all possible values of n in the given range
    for n in range(min_n, max_n + 1):
        # The problem is solvable if and only if n % 7 == 1.
        if n % 7 == 1:
            solvable_n_values.append(n)
            
    # Print the results in increasing order, as requested.
    # The format "each number in the final equation" is interpreted as
    # printing the elements of the solution set.
    print(*solvable_n_values, sep=", ")

find_solvable_hypercube_sizes()