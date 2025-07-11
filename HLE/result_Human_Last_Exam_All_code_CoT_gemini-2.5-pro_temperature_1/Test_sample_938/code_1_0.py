import numpy as np

def solve():
    """
    This function solves the problem by finding the dominant root of the characteristic
    equation and then calculating the required value based on the derived limit.
    """

    # The limit L is given by 2 * ln(lambda), where lambda is the dominant root of
    # the characteristic equation x^4 - 88x - 57 = 0.
    # We will find this root using Newton's method.
    
    # Define the function and its derivative for Newton's method
    def f(x):
        return x**4 - 88*x - 57

    def f_prime(x):
        return 4*x**3 - 88

    # An initial guess for the root. Since x^4 is likely much larger than the other terms,
    # we can approximate x^4 ~ 88x => x^3 ~ 88, so x is around 4.5.
    lambda_root = 4.5
    
    # Iterate Newton's method to get a precise value for the root.
    for _ in range(10):
        lambda_root = lambda_root - f(lambda_root) / f_prime(lambda_root)

    # The limit is L = 2 * ln(lambda)
    limit_val = 2 * np.log(lambda_root)
    
    # We need to compute the integer part of 10^4 * L
    final_value = 10000 * limit_val
    integer_part = int(final_value)

    # Print the steps and the final equation as requested.
    print(f"The characteristic equation of the recurrence is x^4 - 88*x - 57 = 0.")
    print(f"The dominant root, lambda, is found numerically.")
    print(f"lambda approx {lambda_root:.8f}")
    print(f"The limit is L = 2 * ln(lambda) approx {limit_val:.8f}")
    print(f"The value to find is the integer part of 10000 * L.")
    print(f"Final calculation: floor(10000 * 2 * ln({lambda_root:.6f})) = floor({final_value:.6f}) = {integer_part}")

solve()