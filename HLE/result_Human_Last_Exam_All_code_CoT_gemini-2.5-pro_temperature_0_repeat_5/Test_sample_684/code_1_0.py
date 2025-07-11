def find_possible_n_values():
    """
    Finds and prints the values of n for which it is possible to reach a state
    with a single gift on the hypercube.

    The problem is solvable if and only if n is congruent to 1 or 7 modulo 8.
    This code iterates through the given range and checks this condition.
    """
    possible_n = []
    for n in range(7, 56):
        # The condition for solvability is that n mod 8 is 1 or 7.
        if n % 8 == 1 or n % 8 == 7:
            possible_n.append(n)
    
    # The problem asks for the values in increasing order.
    # The loop naturally finds them in increasing order.
    print("The possible values for n are:")
    # The prompt asks to output each number. We will print them separated by commas.
    print(*possible_n, sep=", ")

find_possible_n_values()