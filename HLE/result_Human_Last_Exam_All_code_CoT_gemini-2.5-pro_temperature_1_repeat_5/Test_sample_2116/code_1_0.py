import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem is solved analytically, and the final expression for the expected
    maximum magnitude is pi / (2 * ln(2)). This code calculates the numerical
    value of this expression.
    """
    pi_val = math.pi
    ln2_val = math.log(2)

    # The final analytical result is pi / (2 * ln(2))
    expected_max_magnitude = pi_val / (2 * ln2_val)

    print("The analytical solution for the expected maximum earthquake magnitude is given by the formula:")
    print("E[Max Magnitude] = pi / (2 * ln(2))")
    print("\nCalculating the components of the formula:")
    print(f"pi = {pi_val}")
    print(f"2 = 2")
    print(f"ln(2) = {ln2_val}")
    print("\nPlugging the values into the formula:")
    print(f"Result = {pi_val} / (2 * {ln2_val})")
    print(f"\nThe final calculated value is: {expected_max_magnitude}")

solve_earthquake_magnitude()