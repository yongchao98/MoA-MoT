import math
from fractions import Fraction

def solve_sum():
    """
    Calculates the sum of 1/n^2 for a special set of integers n.
    The sum is known to be equal to zeta(4) * zeta(6) / zeta(2).
    """
    
    # Known values of the Riemann zeta function as rational numbers
    # zeta(2) = pi^2 / 6
    # zeta(4) = pi^4 / 90
    # zeta(6) = pi^6 / 945
    
    # We represent the values without the pi factors for now
    zeta_2_rational = Fraction(1, 6)
    zeta_4_rational = Fraction(1, 90)
    zeta_6_rational = Fraction(1, 945)

    # The power of pi for each zeta value
    zeta_2_pi_power = 2
    zeta_4_pi_power = 4
    zeta_6_pi_power = 6

    # The sum S = zeta(4) * zeta(6) / zeta(2)
    # Calculate the rational part of the result
    result_rational = (zeta_4_rational * zeta_6_rational) / zeta_2_rational

    # Calculate the final power of pi
    result_pi_power = zeta_4_pi_power + zeta_6_pi_power - zeta_2_pi_power

    # Print the equation with the numbers
    print("The sum is evaluated using the identity:")
    print("Sum = zeta(4) * zeta(6) / zeta(2)")
    print("")
    print("Substituting the known values:")
    print(f"Sum = (({zeta_4_rational.numerator}/{zeta_4_rational.denominator}) * pi^{zeta_4_pi_power}) * (({zeta_6_rational.numerator}/{zeta_6_rational.denominator}) * pi^{zeta_6_pi_power}) / (({zeta_2_rational.numerator}/{zeta_2_rational.denominator}) * pi^{zeta_2_pi_power})")
    print("")
    print("Simplifying the expression gives:")
    print(f"Sum = ({result_rational.numerator}/{result_rational.denominator}) * pi^{result_pi_power}")


solve_sum()