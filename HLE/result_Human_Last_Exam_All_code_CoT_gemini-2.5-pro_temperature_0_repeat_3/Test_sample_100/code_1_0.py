import math
from fractions import Fraction

def solve():
    """
    Calculates the coefficients of the final expression for the integral.
    The integral I is expressed as C1*pi^8 + C2*pi^2 + C3*pi + C4.
    """

    # Term IA: integral of p^7 / (e^p - 1)
    # This is Gamma(8) * zeta(8)
    # Gamma(8) = 7!
    # zeta(8) = pi^8 / 9450
    gamma_8 = math.factorial(7)
    # The coefficient of pi^8 is 7! / 9450
    c1 = Fraction(gamma_8, 9450)

    # Term IC: integral of p / (e^p - 1)
    # This is Gamma(2) * zeta(2)
    # Gamma(2) = 1!
    # zeta(2) = pi^2 / 6
    gamma_2 = math.factorial(1)
    # The coefficient of pi^2 from this term is 1! / 6
    c2_part1 = Fraction(gamma_2, 6)

    # Term ID: integral of p*e^-p / (e^p - 1)
    # This is equal to IC - Gamma(2)
    # So it contributes another (pi^2 / 6) and a constant term -1
    c2_part2 = Fraction(gamma_2, 6)
    c4_part1 = -gamma_2

    # Term IB: integral of (e^(p/4) - e^(-p/4)) / (2*(e^p - 1))
    # This integral evaluates to 2 - pi/2
    c3_part1 = Fraction(-1, 2)
    c4_part2 = 2

    # Combine coefficients
    C1 = c1
    C2 = c2_part1 + c2_part2
    C3 = c3_part1
    C4 = c4_part1 + c4_part2

    print("The integral is a sum of four parts, which can be expressed in the form:")
    print("I = C1 * pi^8 + C2 * pi^2 + C3 * pi + C4\n")
    print("The calculated coefficients are:")
    print(f"C1 (coefficient of pi^8): {C1.numerator}/{C1.denominator}")
    print(f"C2 (coefficient of pi^2): {C2.numerator}/{C2.denominator}")
    print(f"C3 (coefficient of pi): {C3.numerator}/{C3.denominator}")
    print(f"C4 (constant term): {C4}\n")

    print("Thus, the final value of the integral is:")
    print(f"({C1.numerator}/{C1.denominator}) * pi^8 + ({C2.numerator}/{C2.denominator}) * pi^2 + ({C3.numerator}/{C3.denominator}) * pi + {C4}")

solve()