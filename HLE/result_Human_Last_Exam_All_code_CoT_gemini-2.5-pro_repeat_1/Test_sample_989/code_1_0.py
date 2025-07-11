import numpy as np
import sympy as sp

def explain_false_statement_E():
    """
    This function explains why statement E is false by providing a counterexample.
    Statement E: Any strictly convex function has a unique global minimizer.

    We will analyze the function f(x) = e^x.
    1. A function is strictly convex if its second derivative is always positive.
    2. A function has a minimum where its first derivative is zero.
    """
    x = sp.Symbol('x')
    f = sp.exp(x)

    # Calculate first and second derivatives
    f_prime = sp.diff(f, x)
    f_double_prime = sp.diff(f_prime, x)

    print(f"Analyzing the function f(x) = {f}")
    print("-" * 40)

    # Check for convexity
    print(f"The second derivative is f''(x) = {f_double_prime}.")
    print("Since e^x is always positive for any real number x, the function f(x) = e^x is strictly convex.")
    print("-" * 40)

    # Check for a minimizer
    print(f"To find a minimum, we set the first derivative to zero: f'(x) = 0.")
    print(f"The first derivative is f'(x) = {f_prime}.")
    print("The equation e^x = 0 has no solution for any real number x.")
    print("Therefore, the function f(x) = e^x has no critical points and thus no minimum value.")
    print("-" * 40)

    print("Conclusion: We have found a strictly convex function (f(x) = e^x) that does not have a global minimizer. This serves as a counterexample, proving that statement E is false.")

if __name__ == '__main__':
    explain_false_statement_E()