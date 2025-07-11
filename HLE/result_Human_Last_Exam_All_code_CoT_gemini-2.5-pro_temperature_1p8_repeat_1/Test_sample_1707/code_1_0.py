import decimal

def solve_for_T():
    """
    Calculates the time T based on the derived solvability condition.

    The solvability condition for the given boundary-value problem simplifies to:
    alpha = T * (sum_{i=1 to inf} x0^i + sum_{i=1 to inf} y0^i)

    This can be written using the formula for a geometric series as:
    T = alpha / (x0 / (1 - x0) + y0 / (1 - y0))

    Given that x0 and y0 are extremely small, the expression (z / (1 - z)) is
    very accurately approximated by z. The calculation is therefore based on
    the highly accurate simplified formula:
    T = alpha / (x0 + y0)
    """

    # Set a high precision for decimal calculations, though for this
    # approximation, it's not strictly necessary. 50 digits is plenty.
    decimal.getcontext().prec = 50

    # Define the given constants as Decimal objects to handle large exponents
    alpha = decimal.Decimal('1e10000')
    x0 = decimal.Decimal('1e-5000000')
    y0 = decimal.Decimal('1e-5000000')

    # The exact equation for T is: T = alpha / (x0/(1-x0) + y0/(1-y0))
    # Due to the extreme values, direct computation of (1-x0) is impractical.
    # We use the highly accurate approximation T = alpha / (x0 + y0).
    T_approx = alpha / (x0 + y0)

    print("The derived equation for T is: T = alpha / (x0/(1-x0) + y0/(1-y0))")
    print("Due to numerical limitations, we use the extremely accurate approximation: T = alpha / (x0 + y0)")
    print("\nSubstituting the given values into the approximated formula:")
    print(f"alpha = {alpha}")
    # Note: In the final formula T = a / (2x_0), since x_0=y_0
    print(f"x0 = {x0}")
    print(f"y0 = {y0}")
    print("\nFinal equation with substituted numbers:")
    print(f"T = {alpha} / ({x0} + {y0})")
    print("\nResult:")
    print(f"T = {T_approx}")


solve_for_T()