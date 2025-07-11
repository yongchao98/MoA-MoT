import sympy

def solve_probability():
    """
    Calculates the probability of reaching a state where 50% of products are good.
    """
    x = sympy.Symbol('x')

    # The problem starts with W_0 = 2 good and B_0 = 1 defective products.
    # The limiting proportion of good products, X_inf, follows a Beta(2, 1) distribution.
    # The PDF of Beta(2, 1) is f(x) = 2x for x in [0, 1].
    pdf = 2 * x

    # Conditional on X_inf = x, the problem is a gambler's ruin problem starting at 1.
    # The probability of ruin (reaching 0) is P_ruin(x).
    # Case 1: x <= 0.5 (fair or unfavorable game for growing the number of good products)
    # P_ruin(x) = 1. The integrand is 1 * 2x.
    integrand_1 = 1 * pdf
    integral_1 = sympy.integrate(integrand_1, (x, 0, sympy.Rational(1, 2)))

    # Case 2: x > 0.5 (favorable game)
    # P_ruin(x) = (1-x)/x. The integrand is ((1-x)/x) * 2x = 2*(1-x).
    prob_ruin_2 = (1 - x) / x
    integrand_2 = prob_ruin_2 * pdf
    integral_2 = sympy.integrate(integrand_2, (x, sympy.Rational(1, 2), 1))

    # The total probability is the sum of the two integrals.
    total_probability = integral_1 + integral_2

    print("The total probability is the sum of two integrals, based on the conditional probability of reaching the target state.")
    print(f"The first integral (from x=0 to 0.5) is: {float(integral_1):.2f}")
    print(f"The second integral (from x=0.5 to 1) is: {float(integral_2):.2f}")
    print("The final equation for the probability is:")
    print(f"{float(integral_1):.2f} + {float(integral_2):.2f} = {float(total_probability):.2f}")
    print(f"\nThe exact probability, and therefore the least upper bound, is {total_probability}.")

solve_probability()
<<<1/2>>>