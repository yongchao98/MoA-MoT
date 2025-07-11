import sympy

def solve_integral():
    """
    This function calculates the definite integral of (sin(2x) - sin(x))^50 from 0 to pi
    using the sympy library for symbolic mathematics.
    """
    # Define the symbolic variable
    x = sympy.Symbol('x')

    # As derived in the explanation, the complex integral simplifies to a more standard form.
    # We showed that the original integral is equivalent to this one:
    integrand = (sympy.sin(2*x) - sympy.sin(x))**50

    # Perform the definite integration from 0 to pi
    integral_result = sympy.integrate(integrand, (x, 0, sympy.pi))

    # The result is of the form (N/D) * pi. We extract the numerator and denominator.
    # The `as_numer_denom()` method separates the expression into a numerator and denominator.
    # We divide by pi first to isolate the rational coefficient.
    rational_part = integral_result / sympy.pi
    numerator, denominator = rational_part.as_numer_denom()

    # Print the components of the final answer
    print("The final integral evaluates to (Numerator * pi) / Denominator")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print(f"Final Answer as a fraction times pi: {numerator}/{denominator} * pi")
    print(f"The full expression is: {integral_result}")

if __name__ == '__main__':
    solve_integral()
