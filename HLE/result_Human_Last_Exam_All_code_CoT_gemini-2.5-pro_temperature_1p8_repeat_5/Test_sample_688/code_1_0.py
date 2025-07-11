import sys
import math
from fractions import Fraction

def calculate_cn_prefactor(n):
    """
    Calculates and prints the prefactor c_n for a given integer n.
    
    The prefactor c_n appears in the virial expansion of the equation of
    state and corresponds to the contribution from the fully f-connected
    Ree-Hoover diagram. The formula is c_n = - (n-1) / n!.

    Args:
        n (int): The number of nodes in the diagram, must be >= 2.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: Input 'n' must be an integer greater than or equal to 2.", file=sys.stderr)
        return

    numerator = n - 1
    
    try:
        denominator = math.factorial(n)
    except ValueError:
        # This handles negative inputs, though we've already checked n < 2
        print(f"Error: Factorial is not defined for n={n}.", file=sys.stderr)
        return

    # To show the fraction in its simplest form
    simplified_fraction = Fraction(numerator, denominator)
    
    # To show the decimal value
    decimal_value = -float(numerator) / float(denominator)

    # Output each number in the final equation as requested.
    print(f"Calculation for n = {n}:")
    print(f"c_{n} = - (n - 1) / n!")
    print(f"c_{n} = - {numerator} / {denominator}")
    
    # If the fraction was simplified, show the simplified form as well.
    if simplified_fraction.denominator != denominator:
         print(f"c_{n} = - {simplified_fraction.numerator} / {simplified_fraction.denominator}")
         
    print(f"c_{n} = {decimal_value}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>", file=sys.stderr)
        print("Description: Calculates the prefactor c_n for the fully f-connected diagram in the virial series.", file=sys.stderr)
        print("Argument <n> must be an integer greater than or equal to 2.", file=sys.stderr)
        sys.exit(1)
    
    try:
        n_value = int(sys.argv[1])
        calculate_cn_prefactor(n_value)
    except ValueError:
        print(f"Error: Invalid input '{sys.argv[1]}'. Please provide a valid integer for n.", file=sys.stderr)
        sys.exit(1)
