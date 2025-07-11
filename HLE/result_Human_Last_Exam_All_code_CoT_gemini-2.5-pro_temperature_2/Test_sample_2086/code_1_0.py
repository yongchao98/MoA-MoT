import math

def calculate_max_eigenvalues():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    The plan is as follows:
    1.  The problem is reduced to finding an optimal partition of the integer n into isomorphism class sizes.
    2.  Analysis of the matrix's eigenvalues shows that the number of eigenvalues greater than 2 is maximized by partitioning n into parts of size 2 and 3.
    3.  This optimization leads to a formula that depends on n modulo 3.
    4.  A concise mathematical expression for this formula is floor((2*n + 3) / 6).
    5.  This script prompts the user for n, calculates the result using this formula,
        and prints the components of the formula and the final integer result.
    """
    try:
        n_str = input("Enter the number of objects n: ")
        n = int(n_str)
        if n <= 0:
            print("Error: The number of objects n must be a positive integer.")
            return

        # The formula derived is floor((2*n + 3) / 6).
        # We will show the components of this formula.
        numerator_coeff = 2
        numerator_add = 3
        denominator = 6

        result = (numerator_coeff * n + numerator_add) // denominator

        print(f"For n = {n}, the formula to calculate the maximum number of eigenvalues greater than 2 is:")
        print(f"result = floor(({numerator_coeff} * n + {numerator_add}) / {denominator})")
        print(f"result = floor(({numerator_coeff * n} + {numerator_add}) / {denominator})")
        print(f"result = floor({numerator_coeff * n + numerator_add} / {denominator})")
        print(f"result = {result}")

    except ValueError:
        print("Error: Invalid input. Please enter a valid integer for n.")

if __name__ == '__main__':
    calculate_max_eigenvalues()