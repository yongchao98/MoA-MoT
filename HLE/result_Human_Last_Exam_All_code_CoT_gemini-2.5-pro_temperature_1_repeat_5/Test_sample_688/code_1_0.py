import math
from fractions import Fraction

def calculate_c_n(n):
    """
    Calculates the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram with n nodes.

    Args:
        n (int): The number of nodes in the diagram (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    print(f"Calculating the prefactor c_n for n = {n}")
    print("-" * 30)

    if n == 2:
        # The case for n=2 is simple as the only biconnected graph is K_2 itself.
        # B_2 is defined directly from the f_12 function.
        # B_2 = -1/2 * integral(f_12), so Lambda_2 = integral(f_12) and c_2 = -1/2.
        c_frac = Fraction(-1, 2)
        print("For n=2, the virial coefficient B_2 is defined as:")
        print("B_2 = -1/2 * integral(f_12)")
        print("The fully f-connected diagram integral is Lambda_2 = integral(f_12).")
        print(f"Therefore, c_2 = {c_frac.numerator}/{c_frac.denominator}")

    else: # n >= 3
        # Use the derived formula: c_n = -(n-1) * (-2)^(n-3) / n!
        n_minus_1 = n - 1
        n_factorial = math.factorial(n)
        exponent = n - 3
        s_n = (-2)**exponent
        
        numerator = -n_minus_1 * s_n
        denominator = n_factorial
        
        c_frac = Fraction(numerator, denominator)

        print("The general formula for n >= 3 is:")
        print("c_n = -(n-1) * (-2)^(n-3) / n!")
        print("\nFor the specific case of n = {}:".format(n))
        print("Numerator = -({} - 1) * (-2)^({} - 3)".format(n, n))
        print("          = -({}) * (-2)^({})".format(n_minus_1, exponent))
        print("          = -({}) * ({}) = {}".format(n_minus_1, s_n, numerator))
        print("\nDenominator = {}! = {}".format(n, n_factorial))
        print("\nResult:")
        print("c_{} = {} / {}".format(n, numerator, denominator))
        print("   = {}/{}".format(c_frac.numerator, c_frac.denominator))

    print("-" * 30)
    print(f"The final value is c_{n} = {float(c_frac)}")


# --- Main execution ---
# Set the desired value of n here
n_value = 5
calculate_c_n(n_value)

<<<
For n=2, c_2 = -1/2. For n>=3, c_n = -(n-1)*(-2)^(n-3)/n!
>>>