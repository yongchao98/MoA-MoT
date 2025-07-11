import math
from fractions import Fraction

def calculate_and_print_c_n(n_values):
    """
    Calculates and prints the prefactor c_n for given values of n.

    The prefactor c_n is given by the formula: c_n = -(n-1) / n!
    
    Args:
        n_values (list[int]): A list of integers n >= 2.
    """
    print("The system-independent prefactor c_n is given by the formula: c_n = -(n-1)/n!")
    print("-" * 40)

    for n in n_values:
        if not isinstance(n, int) or n < 2:
            print(f"Skipping invalid input n={n}. Virial coefficients B_n start from n=2.")
            continue

        numerator = n - 1
        denominator = math.factorial(n)
        
        # Use the Fraction class for an exact fractional representation
        c_n_fraction = Fraction(-numerator, denominator)

        print(f"For n = {n}:")
        # Print the formula with numbers plugged in
        print(f"c_{n} = -({n}-1) / {n}!")
        # Print the intermediate, unsimplified fraction
        print(f"c_{n} = -{numerator} / {denominator}")
        # Print the final simplified fraction
        print(f"c_{n} = {c_n_fraction}")
        print("-" * 20)

if __name__ == '__main__':
    # We will calculate c_n for n from 2 to 5 as an example.
    example_n_values = [2, 3, 4, 5]
    calculate_and_print_c_n(example_n_values)