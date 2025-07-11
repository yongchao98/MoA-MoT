import math

def solve_series():
    """
    Calculates the value of the series sum_{n=2 to inf} (-1)^n * S_n / n
    based on its analytical solution: 1 + (ln(2))^2 / 2 - pi^2 / 12.
    """
    
    # The value of the sum is given by the analytical formula.
    # Let's calculate each term of the formula: 1 + (ln(2))^2 / 2 - pi^2 / 12
    
    term1 = 1.0
    term2 = (math.log(2)**2) / 2
    term3 = (math.pi**2) / 12
    
    # Calculate the final result
    result = term1 + term2 - term3
    
    # As requested, output each number in the final equation
    print("The value of the sum can be calculated using the equation:")
    print(f"{term1} + {term2} - {term3}")
    print("Which evaluates to:")
    print(result)

solve_series()