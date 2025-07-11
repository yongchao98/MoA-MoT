import sys

def solve_fourier_restriction_exponent():
    """
    Calculates the critical exponent for the Fourier restriction problem
    on the moment curve in R^3.

    The problem asks for the largest possible value of p such that no L^p function
    on R^3 can have its Fourier support lying on the moment curve {(t, t^2, t^3)}.
    This value is given by the formula p = n * (n + 1) / 2 for the moment curve in R^n.
    For this problem, n = 3.
    """
    n = 3
    
    # The formula for the critical exponent
    p = n * (n + 1) / 2
    
    print(f"The problem is set in R^n where n = {n}.")
    print("The largest possible value of p is given by the formula: n * (n + 1) / 2")
    print("Plugging in the value of n:")
    # The user asked to show each number in the final equation.
    # We will show the steps of the calculation.
    n_plus_1 = n + 1
    numerator = n * n_plus_1
    
    # We use integer division // as the result is an integer
    p_int = numerator // 2

    print(f"{n} * ({n} + 1) / 2 = {n} * {n_plus_1} / 2 = {numerator} / 2 = {p_int}")

solve_fourier_restriction_exponent()