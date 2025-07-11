def find_possible_n_values():
    """
    This function finds and prints the values of n (between 7 and 55)
    for which it is possible to reach a state with a single gift on the
    5D hypercube.

    The conditions for possibility are:
    1. n must be an odd number.
    2. The remainder of n divided by 7 must be 1 or 6.
    """
    possible_n = []
    for n in range(7, 56):
        # Check if n is odd
        is_odd = (n % 2 != 0)

        # Check if n mod 7 is 1 or 6
        is_mod_7_valid = (n % 7 == 1) or (n % 7 == 6)

        if is_odd and is_mod_7_valid:
            possible_n.append(n)

    # Print the results in increasing order
    for val in possible_n:
        print(val)

find_possible_n_values()
