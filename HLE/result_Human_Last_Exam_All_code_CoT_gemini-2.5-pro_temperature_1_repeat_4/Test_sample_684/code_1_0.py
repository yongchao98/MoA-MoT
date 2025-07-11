def find_possible_n_values():
    """
    This function finds the values of n for which it is possible to leave only one gift.

    The problem is solvable if and only if a condition on n is met. Let n = 7k + r.
    The condition is that the number of 1s in a characteristic vector is exactly one.
    - If k is even, the number of 1s is r. So we need r=1.
    - If k is odd, the number of 1s is 7-r. So we need 7-r=1, which means r=6.

    This script iterates through n from 7 to 55 and checks this condition.
    """
    possible_n = []
    for n in range(7, 56):
        k = n // 7
        r = n % 7

        is_possible = False
        if k % 2 == 0:  # k is even
            if r == 1:
                is_possible = True
        else:  # k is odd
            if r == 6:
                is_possible = True

        if is_possible:
            possible_n.append(n)

    print("The values of n for which it is possible to reach a state with only one gift are:")
    # The problem asks to "output each number in the final equation", which we interpret
    # as printing the list of solutions.
    print(*possible_n, sep=", ")

find_possible_n_values()
