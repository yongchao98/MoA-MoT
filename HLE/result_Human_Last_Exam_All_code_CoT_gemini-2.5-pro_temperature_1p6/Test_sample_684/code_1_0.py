def find_possible_n_values():
    """
    This function finds all integers n in the range [7, 55] for which it is
    possible to leave just one gift on the hypercube.
    
    The problem can be solved using invariants. By constructing a special set of
    weight functions for the positions on the hypercube, one can show that a
    solution is only possible if n is congruent to 1 modulo 7 (n % 7 == 1).
    
    This script iterates through the given range and checks this condition.
    """
    
    start_n = 7
    end_n = 55
    move_length = 7
    
    possible_n = []
    for n in range(start_n, end_n + 1):
        if n % move_length == 1:
            possible_n.append(n)
            
    print("The possible values for n in increasing order are:")
    # The problem asks to output the numbers of the final solution.
    print(*possible_n, sep=', ')

find_possible_n_values()
