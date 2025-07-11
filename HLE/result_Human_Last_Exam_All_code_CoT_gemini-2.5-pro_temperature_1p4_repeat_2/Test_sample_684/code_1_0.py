def find_possible_n_values():
    """
    Finds and prints the values of n (7 <= n <= 55) for which it is possible
    to reach a state with a single gift on a 5D hypercube.

    The derivation shows that this is possible if and only if two conditions are met:
    1. n must be an odd number.
    2. n must not be a multiple of 7.
    
    This script iterates through the given range for n and prints the values
    that satisfy both conditions.
    """
    
    valid_n = []
    for n in range(7, 56):
        # Condition 1: n must be odd.
        is_odd = (n % 2 != 0)
        
        # Condition 2: n must not be a multiple of 7.
        is_not_multiple_of_7 = (n % 7 != 0)
        
        if is_odd and is_not_multiple_of_7:
            valid_n.append(n)
            
    print(*valid_n)

find_possible_n_values()
