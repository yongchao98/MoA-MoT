def find_possible_n_values():
    """
    Finds and prints the integer values of n in the range [7, 55] for which
    the hypercube problem is solvable.
    
    The condition for solvability is that n mod 7 is either 1 or 6.
    """
    
    possible_n = []
    for n in range(7, 56):
        # Check if n satisfies the condition n % 7 == 1 or n % 7 == 6.
        if n % 7 == 1 or n % 7 == 6:
            possible_n.append(n)
            
    # The problem asks for the values in increasing order, which the loop provides.
    # The instruction "output each number in the final equation" is interpreted as
    # printing the list of solutions.
    print(*possible_n, sep=', ')

find_possible_n_values()