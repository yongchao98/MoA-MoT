import math

def solve():
    """
    This function calculates the minimum value of the sum of set sizes
    for a given n, based on the derived formula.
    """
    # We can use any integer n >= 2 as an example. Let's use n = 10.
    n = 10
    
    # The derived formula for the minimum sum is floor(n^2 / 4) + 2.
    
    # Calculate the first term: floor(n^2 / 4)
    term1 = n**2 // 4
    
    # The second term is a constant
    term2 = 2
    
    # Calculate the final result
    min_value = term1 + term2
    
    # Print the breakdown of the calculation as requested.
    print(f"For n = {n}:")
    print(f"The minimum value is given by the formula: floor(n^2 / 4) + 2")
    print(f"Calculation: floor({n}^2 / 4) + {term2} = {term1} + {term2} = {min_value}")

solve()