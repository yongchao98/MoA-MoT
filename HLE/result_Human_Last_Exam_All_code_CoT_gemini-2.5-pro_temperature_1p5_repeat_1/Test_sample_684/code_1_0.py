def find_possible_n_values():
    """
    Finds and prints the values of n for which it is possible to reach
    a state with exactly one gift in the hypercube.

    The problem is solvable if and only if the dimension size n is not divisible by 7.
    This function iterates through n from 7 to 55 and prints the values that
    satisfy this condition.
    """
    
    n_min = 7
    n_max = 55
    
    possible_n_values = []
    for n in range(n_min, n_max + 1):
        # The condition for the puzzle to be solvable is that n is not divisible by 7.
        if n % 7 != 0:
            possible_n_values.append(n)
            
    # Print all the found values in increasing order.
    # The list is already sorted as we iterate from n_min to n_max.
    print(*possible_n_values)

find_possible_n_values()