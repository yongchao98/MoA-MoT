import math

def solve_series():
    """
    Calculates the value of the series sum_{n=2 to inf} (-1)^n * S_n / n
    based on its analytical solution: (ln(2)^2 / 2) - (pi^2 / 12) + 1.
    """
    
    # The value of the series can be found analytically.
    # The final expression for the sum S is:
    # S = (ln(2)^2 / 2) - (pi^2 / 12) + 1

    # Calculate each term of the expression
    ln2 = math.log(2)
    pi = math.pi

    term1 = (ln2**2) / 2
    term2 = (pi**2) / 12
    term3 = 1

    # Calculate the final result
    result = term1 - term2 + term3

    # Print the equation with the numerical values of each term
    print(f"The sum is given by the expression: (ln(2)^2 / 2) + 1 - (pi^2 / 12)")
    print(f"Let's calculate each number in the final equation:")
    print(f"Term 1: (ln(2)^2 / 2) = {term1}")
    print(f"Term 2: 1 = {term3}")
    print(f"Term 3: (pi^2 / 12) = {term2}")
    print(f"\nFinal Equation: {term1} + {term3} - {term2} = {result}")
    print(f"\nThe value of the sum is: {result}")

solve_series()