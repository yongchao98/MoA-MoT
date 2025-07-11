import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
    The limit is determined by the growth rates of two sequences of solutions.
    """
    
    # Numbers from the characteristic equation for the first branch of solutions
    c1 = 13
    d1 = 165
    
    # Numbers from the characteristic equation for the second branch of solutions
    c2 = 5
    d2 = 21

    # The constant factor for each branch's contribution to the limit
    numerator = 2
    
    # The scaling factor from the problem statement
    factor = 10**4
    
    print(f"The final calculation is based on the expression:")
    print(f"{factor} * ( {numerator} / ln( ({c1} + sqrt({d1})) / 2 ) + {numerator} / ln( ({c2} + sqrt({d2})) / 2 ) )")
    print("\nThe numbers in the final equation are:")
    print(f"Factor: {factor}")
    print(f"Numerator: {numerator}")
    print(f"Branch 1 constants: c1={c1}, d1={d1}")
    print(f"Branch 2 constants: c2={c2}, d2={d2}")
    
    # Growth rate for the first sequence
    lambda_val = (c1 + math.sqrt(d1)) / 2
    
    # Growth rate for the second sequence
    mu_val = (c2 + math.sqrt(d2)) / 2
    
    # The limit value
    limit_val = numerator / math.log(lambda_val) + numerator / math.log(mu_val)
    
    # The final result
    result = factor * limit_val
    
    print(f"\nThe value of the limit expression is: {result}")
    print(f"The integer part of the result is: {int(result)}")

solve()
