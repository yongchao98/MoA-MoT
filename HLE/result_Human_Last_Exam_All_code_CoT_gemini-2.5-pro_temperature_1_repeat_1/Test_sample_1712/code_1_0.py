import math

def solve_series():
    """
    Calculates the value of the series sum_{n=2 to infty} (-1)^n * S_n / n.
    The analytical result is (ln(2))^2 / 2 - pi^2 / 12 + 1.
    """
    
    # Calculate the components of the formula
    ln2 = math.log(2)
    pi = math.pi
    
    # Calculate the terms in the final expression
    term1_val = (ln2**2) / 2
    term2_val = (pi**2) / 12
    constant = 1.0
    
    # Calculate the final result
    result = term1_val - term2_val + constant
    
    # Print the final equation with the numerical values of each term
    print("The value of the sum is given by the expression: (ln(2))^2 / 2 - pi^2 / 12 + 1")
    print("The numerical equation is:")
    print(f"{term1_val} - {term2_val} + {constant} = {result}")

solve_series()