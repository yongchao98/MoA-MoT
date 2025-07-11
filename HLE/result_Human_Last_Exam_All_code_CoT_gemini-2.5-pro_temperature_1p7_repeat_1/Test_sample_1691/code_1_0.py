import numpy as np

def solve():
    """
    Calculates the coefficients for the asymptotic approximation of the integral I(epsilon).
    """

    # The integral is I(epsilon) = integral from 0 to 15 of 1/(epsilon + 9x^5 + 5x^6 + 9x^8) dx.
    # For small x, the denominator is dominated by epsilon + b*x^n.
    b = 9.0
    n = 5.0

    # The approximation has the form I(epsilon) ~ C * epsilon^(-a).
    # The exponent 'a' is determined by the power 'n'.
    # a = (n-1)/n
    a = (n - 1.0) / n

    # The coefficient C is given by the formula C = (b**(-1/n) / n) * (pi / sin(pi/n)).
    # We will calculate each component of C.
    pi_val = np.pi
    sin_val = np.sin(pi_val / n)
    b_pow_val = b**(-1.0 / n)
    
    # Calculate C
    C = (b_pow_val / n) * (pi_val / sin_val)

    print("The analytical approximation for I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) approx C * epsilon**(-a)")
    print("\n-------------------------------------------------")

    print("\nStep 1: Determine the exponent 'a'")
    print(f"The dominant term in the denominator near x=0 is f(x) approx b*x**n, where b={b} and n={n}.")
    print(f"The exponent 'a' is calculated as a = (n - 1) / n.")
    print(f"a = ({n} - 1.0) / {n}")
    print(f"a = {a}")

    print("\n-------------------------------------------------")
    print("\nStep 2: Calculate the coefficient C")
    print("The formula for C is: C = (b**(-1/n) / n) * (pi / sin(pi/n))")
    print("\nIntermediate values:")
    print(f"  b = {b}")
    print(f"  n = {n}")
    print(f"  pi = {pi_val}")
    print(f"  sin(pi/n) = sin(pi/{n}) = {sin_val}")
    print(f"  b**(-1/n) = {b}**(-1.0/{n}) = {b_pow_val}")

    print("\nFinal calculation for C:")
    print(f"C = ({b_pow_val} / {n}) * ({pi_val} / {sin_val})")
    print(f"C = {C}")
    
    print("\n-------------------------------------------------")
    print("\nFinal Result:")
    print("The derived analytical formula is:")
    print(f"I(epsilon) approx {C} * epsilon**(-{a})")

solve()