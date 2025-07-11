import math
from fractions import Fraction

def calculate_cn(n_values):
    """
    Calculates the prefactor c_n for a given list of n values.
    The formula is c_n = -(n-1) / n!
    """
    print("Calculating the prefactor c_n = -(n-1) / n! for various n.")
    print("-" * 50)
    for n in n_values:
        if not isinstance(n, int) or n < 2:
            print(f"For n = {n}: Input must be an integer greater than or equal to 2.\n")
            continue

        numerator = -(n - 1)
        denominator = math.factorial(n)
        
        # Use Fraction for an exact result
        cn_frac = Fraction(numerator, denominator)
        
        print(f"For n = {n}:")
        # Outputting each number in the final equation
        print(f"c_{n} = -({n} - 1) / {n}!")
        print(f"c_{n} = {numerator} / {denominator}")
        print(f"c_{n} = {cn_frac}")
        print("-" * 50)

# We can now run the calculation for a few example values of n
example_n = [2, 3, 4, 5, 6]
calculate_cn(example_n)