import math

def propose_equation():
    """
    This function defines the physical parameters and prints the
    final equation for the required velocity v0 with all numbers substituted.
    """
    # Given parameters
    m1 = 1  # kg
    m2 = 2  # kg
    R = 3   # meters
    d = 1   # meters
    g = 10  # m/s^2

    # The general equation is derived from conservation of angular momentum and energy:
    # v0 = (2 / (m2 * d)) * sqrt((m1*R^2 + m2*d^2) * (m1*g*R + m2*g*d))

    print("The proposed equation for the value that v0 must have is constructed by substituting the given values into the derived formula:")

    # Using an f-string to build the equation text.
    # The `^` is used for exponentiation for better readability in the output.
    equation_str = (
        f"v0 = (2 / ({m2} * {d})) * "
        f"sqrt([({m1} * {R}^2 + {m2} * {d}^2) * ({m1} * {g} * {R} + {m2} * {g} * {d})])"
    )

    print(equation_str)

propose_equation()