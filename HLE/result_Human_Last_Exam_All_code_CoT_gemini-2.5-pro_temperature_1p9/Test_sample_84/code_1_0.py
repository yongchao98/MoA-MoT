import math

def solve_for_alpha():
    """
    This function calculates the value of alpha based on the asymptotic
    analysis of the polynomial degree d_n.

    The problem states that for a natural number n, d_n is the smallest degree of a
    real polynomial p_n satisfying:
    1. p_n(i) is in [0, 1] for i in {1, 2, ..., n^2}
    2. p_n(i) is in [2, 3] for i in {n^2+1, ..., n^10}

    The asymptotic growth of d_n is given by d_n = Theta(n^alpha).
    """

    # Step 1: Identify the asymptotic lengths of the two intervals.
    # The first interval's length, L_1, is of the order n^a.
    # I_1 = [1, n^2], so L_1 ~ n^2.
    a = 2

    # The second interval's length, L_2, is of the order n^b.
    # I_2 = [n^2+1, n^10], so L_2 ~ n^10.
    b = 10

    # Step 2: Use the result from Approximation Theory.
    # The minimum degree d_n required to separate two intervals of lengths
    # L_1 ~ n^a and L_2 ~ n^b is given by d_n = Theta(sqrt(L_1 * L_2)).
    #
    # d_n = Theta(sqrt(n^a * n^b))
    #     = Theta(sqrt(n^(a+b)))
    #     = Theta(n^((a+b)/2))
    #
    # We are given d_n = Theta(n^alpha), so alpha = (a+b)/2.

    # Step 3: Calculate alpha and print the result.
    alpha = (a + b) / 2
    
    print("The asymptotic growth rate of d_n is given by d_n = Theta(n^alpha).")
    print("This analysis is based on a theorem for polynomial approximation on two disjoint intervals.")
    print("The first interval has a length proportional to n^a, and the second has a length proportional to n^b.")
    print(f"From the problem description, we have a = {a}.")
    print(f"From the problem description, we have b = {b}.")
    print("\nThe formula for alpha is derived from the geometric mean of the interval lengths:")
    print("d_n = Theta(sqrt(n^a * n^b)) = Theta(n^((a+b)/2))")
    print("Therefore, alpha = (a + b) / 2.")
    print("\nCalculation:")
    print(f"alpha = ({a} + {b}) / 2 = {int(alpha)}")

solve_for_alpha()