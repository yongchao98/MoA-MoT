import math

def solve():
    """
    This function solves the problem by:
    1. Defining the characteristic polynomial f(x) and its derivative f'(x).
    2. Using Newton's method to find the dominant root lambda.
    3. Calculating the final value 10^4 * ln(lambda).
    4. Printing the integer part of the result.
    """
    
    # The limit L is ln(lambda), where lambda is the largest root of x^4 - 88x - 57 = 0.
    f = lambda x: x**4 - 88*x - 57
    f_prime = lambda x: 4*x**3 - 88

    # Use Newton's method to find the root.
    # From f(4) < 0 and f(5) > 0, we can start with an initial guess of 4.5.
    lambda_val = 4.5
    for _ in range(10):  # 10 iterations are more than enough for high precision
        lambda_val = lambda_val - f(lambda_val) / f_prime(lambda_val)

    # The problem asks for the integer part of 10^4 * L
    limit_val = math.log(lambda_val)
    result = 10**4 * limit_val
    
    # Output the final calculation and its integer part
    print(f"The dominant root lambda is approximately: {lambda_val}")
    print(f"The limit is ln(lambda), which is approximately: {limit_val}")
    print(f"The value to find is 10^4 * ln(lambda):")
    equation_str = f"10000 * {limit_val:.6f} = {result:.6f}"
    print(equation_str)
    print(f"The integer part of the result is: {math.floor(result)}")

solve()