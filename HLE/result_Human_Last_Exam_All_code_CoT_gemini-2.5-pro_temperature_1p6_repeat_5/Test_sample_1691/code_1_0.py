import math

def approximate_integral():
    """
    Calculates the constants for the analytical approximation of the integral I(epsilon)
    for the small epsilon regime.
    """
    
    # The integral is of the form I(epsilon) = integral[1 / (epsilon + f(x))] dx.
    # For small epsilon, the behavior is dominated by the region where x is small.
    # The function f(x) = 9.0*x**5.0 + 5.0*x**6.0 + 9.0*x**8.0 is approximated
    # by its dominant term for x->0, which is c*x**k.
    c = 9.0
    k = 5.0
    
    # The integral approximation for small epsilon is I(epsilon) ≈ C * epsilon**(-p).
    
    # The power 'p' is determined by the exponent 'k'.
    # p = (k-1)/k
    p = (k - 1.0) / k
    
    # The coefficient 'C' is found from the analytical solution of the simplified integral:
    # I ≈ integral from 0 to infinity of 1/(epsilon + c*x**k) dx
    # This evaluates to C = (π / (k * sin(π/k))) / c**(1/k)
    coefficient_C = math.pi / (k * math.pow(c, 1.0/k) * math.sin(math.pi/k))

    print("The integral is I(e) = integral_0^15 [1 / (e + 9x^5 + 5x^6 + 9x^8)] dx.")
    print("For small epsilon (e), the integral is approximated by I(e) ≈ C * e**(-p).")
    print("\n-----------------------------------------------------")
    print("The final derived analytical formula is:")
    print(f"I(epsilon) ≈ {coefficient_C:.5f} * epsilon**(-{p:.1f})")
    print("-----------------------------------------------------\n")
    print("The individual numbers in this equation are:")
    print(f"Coefficient C = {coefficient_C:.5f}")
    print(f"Exponent p = {p:.1f}")

# Execute the function to print the results
approximate_integral()