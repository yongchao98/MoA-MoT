import math

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # Grid size
    n = 5

    # The problem provides the positions of 'o's and 'x's. We represent these
    # as permutations. pi_O[i] is the row of the 'o' in column i.
    # Positions of 'o's: (1,1), (2,2), (3,3), (4,4), (5,5)
    # We use 1-based indexing for our lists by adding a dummy 0 at the start.
    pi_O = [0, 1, 2, 3, 4, 5]

    # Positions of 'x's: (1,4), (2,5), (3,1), (4,2), (5,3)
    pi_X = [0, 4, 5, 1, 2, 3]

    print(f"The grid size is n = {n}.")
    print(f"The permutation for O markers is: {pi_O[1:]}")
    print(f"The permutation for X markers is: {pi_X[1:]}")

    writhe = 0
    # The writhe is calculated by the formula:
    # w = sum_{1 <= i < j <= n} sgn((pi_O[i] - pi_O[j]) * (pi_X[i] - pi_X[j]))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            # sgn(x) is 1 if x > 0, -1 if x < 0, and 0 if x == 0.
            term_sign = (pi_O[i] - pi_O[j]) * (pi_X[i] - pi_X[j])
            if term_sign > 0:
                writhe += 1
            elif term_sign < 0:
                writhe -= 1

    print(f"\nThe writhe (w) of the diagram is calculated as {writhe}.")

    # The maximal Thurston-Bennequin number (tb) is given by tb = w - n.
    tb = writhe - n

    print("\nThe maximal Thurston-Bennequin number (tb) is found using the formula: tb = w - n")
    # We print each number in the final equation.
    print(f"tb = {writhe} - {n}")
    print(f"The result is tb = {tb}")

solve()