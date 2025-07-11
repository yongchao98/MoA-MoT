import math

def approximate_integral():
    """
    This function develops and prints an analytical formula that approximates the integral
    I(epsilon) for a small epsilon regime.
    """
    
    # Parameters from the dominant term of the denominator, f(x) approx a*x^p
    a = 9.0
    p = 5.0
    
    # The integral is of the form I(epsilon) = integral from 0 to L of 1 / (epsilon + f(x)) dx.
    # For small epsilon, this is approximated by I(epsilon) ~ Integral(0, inf) dx / (epsilon + a*x^p).
    # The solution to this approximate integral is C * epsilon^k.
    
    # 1. Calculate the exponent k
    k = 1.0 / p - 1.0
    
    # 2. Calculate the coefficient C
    # C = (pi / (p * sin(pi/p))) * a^(-1/p)
    C = (math.pi / (p * math.sin(math.pi / p))) * (a**(-1.0 / p))
    
    # 3. Print the derivation and the final formula
    print("The analytical approximation for the integral I(\u03B5) in the small \u03B5 regime is of the form:")
    print("I(\u03B5) \u2248 C * \u03B5^k")
    print("\nThe constants are derived from the dominant term of the denominator, 9.0 * x^5.0, where a=9.0 and p=5.0.")
    
    print("\nThe exponent k is calculated as:")
    # Using ':.1f' to represent the floating point numbers from the problem statement
    print(f"k = 1/p - 1 = 1/{p:.1f} - 1 = {k:.1f}")
    
    print("\nThe coefficient C is calculated using the formula:")
    print(f"C = (\u03C0 / (p * sin(\u03C0/p))) * a^(-1/p)")
    # Showing the values plugged into the formula
    print(f"C = (\u03C0 / ({p:.1f} * sin(\u03C0/{p:.1f}))) * {a:.1f}^(-1/{p:.1f})")
    print(f"C \u2248 {C:.17f}")
    
    print("\n--------------------------------------------------")
    print("The final approximate analytical formula is:")
    # Using unicode for epsilon and approx symbol for better readability
    print(f"I(\u03B5) \u2248 {C:.5f} * \u03B5^({k:.1f})")
    print("--------------------------------------------------")

approximate_integral()