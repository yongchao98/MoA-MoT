import math

def solve():
    """
    Finds the integer part of 10^4 * lim_{n->inf} (ln(s_n))/n.
    The limit is ln(lambda), where lambda is the dominant root of x^4 - 88x - 57 = 0.
    """

    # Define the function and its derivative for Newton's method
    def f(x):
        return x**4 - 88*x - 57

    def f_prime(x):
        return 4*x**3 - 88

    # Use Newton's method to find the root lambda.
    # From analysis, we know there's a root between 4 and 5.
    # Let's start with an initial guess.
    lambda_val = 4.5
    # Iterate to get a precise value for lambda
    for _ in range(10):
        lambda_val = lambda_val - f(lambda_val) / f_prime(lambda_val)

    # The factor from the problem statement
    factor = 10000

    # Calculate the natural logarithm of lambda
    log_lambda = math.log(lambda_val)

    # Calculate the final result
    result = factor * log_lambda

    # Print the details of the final calculation
    print(f"The equation for the recurrence relation is a_n = 88*a_(n-3) + 57*a_(n-4)")
    print(f"The dominant root of the characteristic equation x^4 - 88*x - 57 = 0 is found to be lambda.")
    print(f"Calculated lambda: {lambda_val}")
    print(f"The limit expression is {factor} * ln(lambda).")
    print(f"{factor} * ln({lambda_val}) = {result}")
    
    # Get the integer part of the result
    integer_part = int(result)
    
    print(f"The integer part of the result is {integer_part}.")


solve()
<<<15363>>>