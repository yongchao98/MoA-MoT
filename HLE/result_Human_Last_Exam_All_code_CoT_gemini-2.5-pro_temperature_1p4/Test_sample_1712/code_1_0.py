import math

def solve_sum():
    """
    Calculates the value of the infinite sum S = sum_{n=2 to inf} (-1)^n * S_n / n
    where S_n is the n-th harmonic number.

    The sum can be analytically solved to the expression: 1 - (pi^2 / 12) + (ln(2)^2 / 2).
    This function computes the value of this expression.
    """
    
    # Define constants from the math module
    pi = math.pi
    ln2 = math.log(2)

    # Calculate the three terms in the final expression
    term1 = 1.0
    term2 = (pi**2) / 12
    term3 = (ln2**2) / 2

    # Calculate the final result
    result = term1 - term2 + term3

    # Print the derivation steps as requested
    print("The value of the sum can be shown to be equal to the expression:")
    print("S = 1 - (pi^2 / 12) + (ln(2)^2 / 2)")
    print("\nLet's break down the calculation:")
    
    print(f"The first term is: {term1}")
    print(f"The second term, pi^2 / 12, is: {pi**2} / 12 = {term2}")
    print(f"The third term, (ln(2))^2 / 2, is: {ln2**2} / 2 = {term3}")

    print("\nSubstituting these values into the equation:")
    print(f"S = {term1} - {term2} + {term3}")
    
    print("\nThe final numerical value is:")
    print(result)

solve_sum()