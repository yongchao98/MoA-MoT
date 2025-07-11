import math

def solve():
    """
    This function calculates the integer part of the limit expression
    10^4 * lim_{N->inf} F(N)/ln(N).
    The derivation shows the limit is 2 / ln(alpha), where alpha is a root of x^2-13x+1=0.
    """
    
    # Constants from the derived formula for the limit
    c1 = 10000
    c2 = 2
    c3 = 13  # from the recurrence relation u_{n+1} = 13*u_n - u_{n-1} - 5
    c4 = 165 # from the discriminant of the characteristic equation r^2 - 13r + 1 = 0
    c5 = 2
    
    # Calculate alpha, the dominant root of the characteristic equation
    alpha = (c3 + math.sqrt(c4)) / c5
    
    # Calculate the limit value
    limit_val = c2 / math.log(alpha)
    
    # Calculate the final result
    result = c1 * limit_val
    
    # Print the equation and the values of its components as requested
    print("The final expression is: c1 * (c2 / ln((c3 + sqrt(c4)) / c5))")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    
    # Print the integer part of the result
    print(f"The integer part of the result is: {int(result)}")

solve()
