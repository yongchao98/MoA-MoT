import math
from fractions import Fraction

def calculate_prefactor_cn(n):
    """
    This function calculates the system-independent prefactor c_n for a given integer n >= 2.

    The prefactor c_n is the coefficient of the fully f-connected diagram in the
    virial expansion for the n-th virial coefficient B_n.

    The formula derived is: c_n = -(n-1) / n!
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    print(f"Calculating the prefactor c_n for n = {n}.")
    print(f"The formula is: c_n = -(n - 1) / n!")
    print("-" * 30)

    # Calculate the components of the formula
    numerator_val = n - 1
    denominator_val = math.factorial(n)

    # Use the fractions module for an exact, simplified representation
    result_fraction = Fraction(-numerator_val, denominator_val)

    # Output the final equation with each number explicitly shown
    print("The final equation is:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"   = -{numerator_val} / {denominator_val}")
    print(f"   = {result_fraction.numerator} / {result_fraction.denominator}")
    print("-" * 30)


# --- Main execution ---
# Let's calculate the prefactor for a specific value, n = 4, as a demonstration.
if __name__ == "__main__":
    n_value = 4
    calculate_prefactor_cn(n_value)
    # The prefactor for n=4 is c_4 = -(4-1)/4! = -3/24 = -1/8
    final_answer = Fraction(-(n_value-1), math.factorial(n_value))
    print(f"<<<c_{n_value} = {final_answer.numerator}/{final_answer.denominator}>>>")