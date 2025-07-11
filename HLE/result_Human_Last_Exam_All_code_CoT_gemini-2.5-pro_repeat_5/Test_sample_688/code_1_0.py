import math
from fractions import Fraction

def get_cn(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    Ree-Hoover diagram with n nodes.

    Args:
        n (int): The number of nodes (must be >= 2).

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")

    if n == 2:
        # The case n=2 is an exception to the general formula.
        # As derived from first principles: c_2 = w_K2 * (-1/2!) = 1 * (-1/2) = -1/2
        return -0.5
    else:
        # For n >= 3, the general formula is c_n = (-1)^n / n.
        return ((-1)**n) / n

def main():
    """
    Main function to explain and print the prefactors c_n.
    """
    print("The system-independent prefactor c_n for the fully f-connected Ree-Hoover diagram is determined as follows:")
    print("The prefactor depends on the number of nodes, n.")
    print("\nFor n = 2, the prefactor is a special case:")
    n = 2
    c_n = get_cn(n)
    fraction_c_n = Fraction(c_n).limit_denominator()
    print(f"c_{n} = {fraction_c_n}")

    print("\nFor n >= 3, the prefactor follows the general formula c_n = (-1)^n / n:")
    for n in range(3, 8):
        c_n = get_cn(n)
        fraction_c_n = Fraction(c_n).limit_denominator()
        # The prompt requires outputting each number in the final equation.
        # We will show the calculation for clarity.
        sign_str = "-" if n % 2 != 0 else ""
        print(f"c_{n} = (-1)^{n} / {n} = {sign_str}1/{n} = {fraction_c_n}")

if __name__ == "__main__":
    main()
