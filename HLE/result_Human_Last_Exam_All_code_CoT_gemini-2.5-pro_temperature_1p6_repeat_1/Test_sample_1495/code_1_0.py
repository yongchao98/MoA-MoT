import numpy as np

def solve_probability():
    """
    This function calculates the probability that a randomly chosen point P
    in triangle ABC lies within the inner triangle XYZ, which is formed by
    the intersections of cevians AD, BE, and CF.
    Points D, E, and F are chosen uniformly at random on the sides BC, AC,
    and AB, respectively.

    The probability is the expected value of the ratio Area(XYZ)/Area(ABC).
    We use a Monte Carlo simulation to estimate this value.
    """

    # Set the number of random samples for the simulation.
    # A larger number gives a more accurate result.
    N = 10**7

    # Generate N random triples (u, v, w), where u, v, w are drawn from a
    # uniform distribution over [0, 1].
    # u represents BD/BC, v represents CE/CA, w represents AF/AB.
    rand_vars = np.random.rand(N, 3)
    u, v, w = rand_vars[:, 0], rand_vars[:, 1], rand_vars[:, 2]

    # Calculate the area ratio based on Routh's Theorem.
    # The formula for the ratio Area(XYZ)/Area(ABC) is:
    # ((1-u)(1-v)(1-w) - uvw)^2 / ((1-v+uv)(1-w+vw)(1-u+wu))

    # Calculate the numerator of the ratio formula
    term1_num = (1 - u) * (1 - v) * (1 - w)
    term2_num = u * v * w
    numerator = (term1_num - term2_num)**2

    # Calculate the denominator of the ratio formula
    term1_den = 1 - v + u * v
    term2_den = 1 - w + v * w
    term3_den = 1 - u + w * u
    denominator = term1_den * term2_den * term3_den

    # The area ratio for each random configuration.
    # We must handle the cases where the denominator is zero. In these
    # specific cases, the numerator is also zero, and the limit is 0.
    # The np.divide function with a 'where' clause or post-hoc replacement
    # of NaN with 0 can handle this.
    ratios = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=(denominator != 0))

    # The probability is the expected value of the ratios, which we approximate
    # with the mean of our sampled ratios.
    probability = np.mean(ratios)

    print(f"The estimated probability is: {probability}")

solve_probability()