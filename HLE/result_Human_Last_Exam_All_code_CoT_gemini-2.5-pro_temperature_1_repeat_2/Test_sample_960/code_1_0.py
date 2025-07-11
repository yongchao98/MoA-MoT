import sympy

def solve_probability():
    """
    Calculates the probability of reaching a state with 50% good products.
    The method is based on integrating the conditional probability (from a gambler's ruin problem)
    over the limit distribution of the proportion of good products (a Beta distribution).
    """

    y = sympy.Symbol('y')

    # PDF of the Beta(2, 1) distribution
    pdf = 2 * y

    # Conditional probability of stopping
    # g(y) = 1 if y <= 1/2, (1-y)/y if y > 1/2
    # The integral is split into two parts based on this condition.

    # Part 1: Integral from 0 to 1/2
    integrand1 = 1 * pdf
    integral1 = sympy.integrate(integrand1, (y, 0, 1/2))

    # Part 2: Integral from 1/2 to 1
    integrand2 = ((1 - y) / y) * pdf
    # Note: the 'y' in the denominator of the conditional probability cancels out with a 'y' from the pdf
    # So the integrand becomes 2*(1-y)
    simplified_integrand2 = 2 * (1 - y)
    integral2 = sympy.integrate(simplified_integrand2, (y, 1/2, 1))

    # Total probability
    total_prob = integral1 + integral2

    # Print the explanation and calculation steps
    print("The probability P is the sum of two integrals:")
    print(f"P = \u222B[0, 1/2] (1 * ({pdf})) dy + \u222B[1/2, 1] (((1-y)/y) * ({pdf})) dy")
    print(f"P = \u222B[0, 1/2] {integrand1} dy + \u222B[1/2, 1] {simplified_integrand2} dy")
    print("\nCalculating the first integral:")
    print(f"\u222B {integrand1} dy = {sympy.integrate(integrand1, y)}")
    print(f"Value = [{sympy.integrate(integrand1, y)}]_0^(1/2) = ({sympy.sympify(1/2)**2}) - (0) = {integral1}")

    print("\nCalculating the second integral:")
    print(f"\u222B {simplified_integrand2} dy = {sympy.integrate(simplified_integrand2, y)}")
    print(f"Value = [{sympy.integrate(simplified_integrand2, y)}]_ (1/2)^1 = (2*1 - 1**2) - (2*(1/2) - (1/2)**2) = 1 - (1 - 1/4) = {integral2}")

    print("\nThe total probability is the sum of the two parts:")
    print(f"P = {integral1} + {integral2} = {total_prob}")
    print(f"\nThe upper bound for the probability is {total_prob}.")

solve_probability()