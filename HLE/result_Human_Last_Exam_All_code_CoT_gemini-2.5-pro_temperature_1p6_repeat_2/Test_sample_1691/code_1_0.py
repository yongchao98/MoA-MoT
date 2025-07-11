import numpy as np

def solve_integral_approximation():
    """
    This function calculates and prints the analytical approximation for the integral
    I(epsilon) for small epsilon.
    """
    
    # Define the constants from the denominator polynomial g(x) = c5*x^5 + c6*x^6 + ...
    c5 = 9.0
    c6 = 5.0
    
    # The dominant power of x in g(x) for small x
    n = 5.0
    
    # The powers of epsilon in the resulting asymptotic expansion are derived from scaling analysis.
    # The leading term power is p1 = (n-1)/n
    # The first correction term power is p2 = (n-2)/n
    p1 = (n - 1.0) / n
    p2 = (n - 2.0) / n
    
    # --- Calculate the coefficient K1 for the leading term ---
    # K1 = c5**(-1/n) * Integral(1/(1+u**n) du) from 0 to infinity
    # The value of the definite integral is (pi/n) / sin(pi/n)
    integral1_val = (np.pi / n) / np.sin(np.pi / n)
    K1 = c5**(-1.0/n) * integral1_val
    
    # --- Calculate the coefficient K2 for the first correction term ---
    # K2 = -c6 * c5**(-(n+2)/n) * Integral(u**(n+1)/(1+u**n)**2 du)
    # The value of this integral can be shown to be (2/n^2) * pi / sin(2*pi/n)
    integral2_val = (2.0 / n**2) * np.pi / np.sin(2.0 * np.pi / n)
    K2 = -c6 * c5**(-(n + 2.0) / n) * integral2_val
    
    # Print the derived analytical formula with the calculated coefficients.
    # The final equation is I(epsilon) = K1 * epsilon**(-p1) + K2 * epsilon**(-p2)
    print("An analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) = {K1:.6f} * epsilon**(-{p1:.2f}) + ({K2:.6f}) * epsilon**(-{p2:.2f})")
    
    # To make it cleaner for display:
    print("\nOr more cleanly:")
    print(f"I(epsilon) = {K1:.6f} * epsilon**(-{p1:.2f}) - {abs(K2):.6f} * epsilon**(-{p2:.2f})")
    
    print("\nEach number in the final equation is:")
    print(f"Coefficient K1 = {K1}")
    print(f"Exponent p1 = {p1}")
    print(f"Coefficient K2 = {K2}")
    print(f"Exponent p2 = {p2}")


solve_integral_approximation()