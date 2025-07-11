def solve_alpha():
    """
    Calculates the value of alpha based on the analysis of the polynomial degree growth.

    The problem is to find the asymptotic growth rate d_n = Theta(n^alpha) for the smallest
    degree d_n of a real polynomial p_n satisfying certain conditions on two sets of points.

    1. The problem is equivalent to finding the minimal degree of a polynomial that separates
       the continuous intervals I_1 = [1, n^2] and I_2 = [n^2+1, n^{10}].

    2. From approximation theory, the degree d_n required to separate two intervals
       is proportional to sqrt(L1 * L2) / g, where L1 and L2 are the lengths of the
       intervals and g is the gap between them.

    3. We determine the asymptotic orders of these lengths in terms of n.
       - Length of I_1 = n^2 - 1 = Theta(n^2). So, beta_1 = 2.
       - Length of I_2 = n^{10} - n^2 - 1 = Theta(n^{10}). So, beta_2 = 10.
       - Gap g = (n^2+1) - n^2 = 1 = Theta(n^0). So, gamma = 0.

    4. The asymptotic order of the degree d_n is therefore:
       d_n = Theta(sqrt(n^beta_1 * n^beta_2) / n^gamma)
           = Theta(sqrt(n^(beta_1 + beta_2)) / n^gamma)
           = Theta(n^((beta_1 + beta_2)/2 - gamma))

    5. We are given d_n = Theta(n^alpha), so alpha = (beta_1 + beta_2)/2 - gamma.
    """

    # Exponent for the length of the first interval, [1, n^2]
    beta_1 = 2
    # Exponent for the length of the second interval, [n^2+1, n^{10}]
    beta_2 = 10
    # Exponent for the length of the gap between the intervals
    gamma = 0

    # Calculate alpha using the formula derived
    alpha = (beta_1 + beta_2) / 2 - gamma

    print(f"The exponent for the length of the first set of points is beta_1 = {beta_1}")
    print(f"The exponent for the length of the second set of points is beta_2 = {beta_2}")
    print(f"The exponent for the gap between the sets is gamma = {gamma}")
    print(f"The formula for alpha is (beta_1 + beta_2) / 2 - gamma")
    print(f"alpha = ({beta_1} + {beta_2}) / 2 - {gamma} = {int(alpha)}")

solve_alpha()