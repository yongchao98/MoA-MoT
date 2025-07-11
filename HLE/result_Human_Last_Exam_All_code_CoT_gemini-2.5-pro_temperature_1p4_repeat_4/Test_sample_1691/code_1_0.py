import numpy as np

def solve_integral_approximation():
    """
    This function calculates the analytical approximation for the integral I(e)
    for small e and prints the derivation and result.
    """
    # Define the parameters from the integral's denominator
    # g(x) = c5*x^n + c6*x^d + c8*x^e
    c5 = 9.0
    n = 5.0
    
    # The approximation for small epsilon is of the form I(e) ≈ C * e^p.
    # The power p is determined by the dominant term at x=0.
    # p = 1/n - 1
    p = 1.0 / n - 1.0

    # The coefficient C is derived from the analytical evaluation of the integral.
    # C = (integral from 0 to inf of 1/(1+u^n) du) / c5^(1/n)
    # The standard integral evaluates to (pi/n)/sin(pi/n).
    integral_val = (np.pi / n) / np.sin(np.pi / n)
    C = integral_val / (c5**(1.0 / n))

    # Output the final formula with the numbers from the equation.
    print("The analytical approximation for I(e) is derived from the dominant behavior near x=0.")
    print("The formula has the form: I(e) ≈ C * e^p")
    print("\nThe equation for the power 'p' is:")
    print(f"p = 1/{int(n)} - 1 = {p}")

    print("\nThe equation for the coefficient 'C' is:")
    print(f"C = (pi/{int(n)}) / (sin(pi/{int(n)}) * ({c5})^(1/{int(n)}))")
    
    print("\nPlugging in the values, the final equation is:")
    # Using f-string formatting to display the numbers in the formula
    print(f"I(e) ≈ ( (pi/{int(n)}) / (sin(pi/{int(n)}) * {c5}**(1/{int(n)})) ) * e^(1/{int(n)} - 1)")
    
    # Print the evaluated numerical formula
    print("\nThe evaluated numerical approximation is:")
    print(f"I(e) ≈ {C:.5f} * e^({p})")

solve_integral_approximation()