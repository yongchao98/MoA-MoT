def find_possible_n_values():
    """
    Finds and prints the values of n for which it is possible to reach the state with only one gift.
    The conditions are:
    1. 7 <= n <= 55
    2. n must be odd.
    3. n must be congruent to 1 modulo 7 (n % 7 == 1).
    """
    possible_n = []
    for n in range(7, 56):
        # n must be odd
        is_odd = (n % 2 != 0)
        # n must be 1 mod 7
        is_one_mod_seven = (n % 7 == 1)

        if is_odd and is_one_mod_seven:
            possible_n.append(n)

    for n_val in possible_n:
        print(n_val)

find_possible_n_values()