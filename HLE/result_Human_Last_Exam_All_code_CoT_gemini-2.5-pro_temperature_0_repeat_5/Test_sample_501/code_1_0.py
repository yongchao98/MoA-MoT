import math

def print_force_law():
    """
    This function prints the derived force law for a thermally isolated polymer chain.
    
    The force F is a function of the extension x, the number of segments n,
    the segment length l, and the kinetic energy at zero extension E(0).
    """
    
    # The derived formula is:
    # F(x) = - (2 * E(0) * x / (n**2 * l**2)) * exp(x**2 / (n**2 * l**2))
    
    print("The derived force law between the polymer ends is:")
    print("F(x) = - (2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))")
    print("\nWhere:")
    print("  F(x) is the restoring force at extension x.")
    print("  E(0) is the kinetic energy of the polymer at zero extension.")
    print("  x is the separation of the polymer ends.")
    print("  n is the number of segments (assuming n is large).")
    print("  l is the length of each segment.")
    print("  exp() is the exponential function.")

    print("\nTo meet the output requirements, here are the numerical constants in the final equation:")
    # The equation can be written as F(x) = - (A * E(0) * x / (n**B * l**C)) * exp(x**D / (n**E * l**F))
    A = 2
    B = 2
    C = 2
    D = 2
    E = 2
    F = 2
    
    print(f"The coefficient of the pre-exponential factor is: {A}")
    print(f"The exponent of 'n' in the pre-exponential denominator is: {B}")
    print(f"The exponent of 'l' in the pre-exponential denominator is: {C}")
    print(f"The exponent of 'x' in the exponential term is: {D}")
    print(f"The exponent of 'n' in the exponential denominator is: {E}")
    print(f"The exponent of 'l' in the exponential denominator is: {F}")

# Execute the function to print the result
print_force_law()