import math

def solve_trajectory():
    """
    Calculates the final position of the particle.
    """
    # The final position x(t) is given by the equation:
    # x(t) = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    # Here are the numbers used in this equation:
    n1 = 18.0
    n2 = 6.0
    n3 = 3.0
    power = 1.0 / 3.0

    # Calculate the value of sqrt(3)
    sqrt_3 = math.sqrt(n3)

    # Calculate the terms inside the cube roots
    radicand1 = n1 - n2 * sqrt_3
    radicand2 = n1 + n2 * sqrt_3

    # Calculate the cube roots
    term1 = radicand1**power
    term2 = radicand2**power

    # The final position is the sum of these two terms
    final_position = term1 + term2

    print("The final equation for the position x(t) is:")
    print(f"x(t) = ({n1} - {n2}*sqrt({n3}))**({power:.3f}) + ({n1} + {n2}*sqrt({n3}))**({power:.3f})")
    print("\nThe numerical value of the final position is:")
    print(final_position)

solve_trajectory()