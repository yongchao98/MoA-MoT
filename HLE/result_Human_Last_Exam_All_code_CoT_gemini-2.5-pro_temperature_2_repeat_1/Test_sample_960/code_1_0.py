import sympy

def solve():
    """
    This script calculates the least upper bound for the probability that the production
    process reaches a point where exactly 50% of the products are good.

    The method is as follows:
    1. The process is a Polya's Urn model, which is equivalent to a Bernoulli process
       with a success probability 'p' drawn from a Beta distribution.
    2. The initial state (W=2, B=1) implies a Beta(2, 1) distribution for 'p'.
       The PDF is f(p) = 2*p.
    3. The condition for stopping (W_t = B_t) is equivalent to a 1D random walk
       (representing new products) hitting the level -1.
    4. The probability of this random walk hitting -1, given 'p', is:
       - 1, if p <= 1/2
       - (1-p)/p, if p > 1/2
    5. The total probability is the integral of this hitting probability weighted by the
       PDF of p, from p=0 to 1.
    """

    # Define the symbolic variable for probability
    p = sympy.Symbol('p')

    # PDF of the Beta(2, 1) distribution
    pdf = 2 * p

    # Hitting probability function pi(p)
    # This is split into two parts for the integration.

    # Integral for the case p <= 1/2, where hitting probability is 1
    # We integrate from p=0 to p=1/2
    integral_part1 = sympy.integrate(pdf * 1, (p, 0, sympy.Rational(1, 2)))

    # Integral for the case p > 1/2, where hitting probability is (1-p)/p
    # We integrate from p=1/2 to p=1
    hitting_prob_part2 = (1 - p) / p
    integral_part2 = sympy.integrate(pdf * hitting_prob_part2, (p, sympy.Rational(1, 2), 1))

    # The total probability is the sum of these two parts.
    # This exact probability is the least upper bound.
    total_probability = integral_part1 + integral_part2

    # Print the step-by-step calculation and the final result
    print("The total probability is the sum of two integrals based on the value of p.")
    print(f"The first integral, from p=0 to 1/2, evaluates to: {integral_part1}")
    print(f"The second integral, from p=1/2 to 1, evaluates to: {integral_part2}")
    print("The final equation is:")
    print(f"P(T < \u221e) = {integral_part1} + {integral_part2} = {total_probability}")
    print(f"Thus, the least upper bound for the probability is {total_probability}.")

solve()
<<<0.5>>>