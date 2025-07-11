import sympy

def solve_product_inequality():
    """
    This function states the result for the product inequality in SL_2(R).
    The problem is to find the largest K such that for any compact subset X
    of SL_2(R), mu(X^3) >= K * mu(X), where mu is the Haar measure.

    This is a known result from advanced mathematics (additive combinatorics
    on Lie groups). The value of K is not found by a simple computation but
    by deep structural theorems about the group SL_2(R).

    The largest possible value for K is 4.
    """
    K = 4
    
    # Using sympy to represent the mathematical objects for clarity
    mu = sympy.Function('mu')
    X = sympy.Symbol('X')
    X_cubed = sympy.Symbol('X^3')
    
    # The inequality is mu(X^3) >= K * mu(X)
    inequality = sympy.Ge(mu(X_cubed), K * mu(X))
    
    # We want to print the equation with the variables and the determined value of K.
    # The variables are mu(X^3), K, and mu(X)
    print(f"The inequality is: {mu(X_cubed)} >= {K}*{mu(X)}")
    print(f"The largest possible value of K is {K}.")

solve_product_inequality()