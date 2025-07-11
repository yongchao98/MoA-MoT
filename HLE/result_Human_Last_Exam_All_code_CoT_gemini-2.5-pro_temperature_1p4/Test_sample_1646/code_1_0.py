import numpy as np
from scipy.special import gamma

def f(z):
    """
    The derived function f(z) = sqrt(pi) / Gamma(z).
    """
    return np.sqrt(np.pi) / gamma(z)

def verify_functional_equation(z):
    """
    Verifies the functional equation for a given complex number z.
    The equation is f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2).
    """
    print(f"--- Verifying for z = {z} ---")
    
    # Calculate Left-Hand Side (LHS)
    lhs = f(z)
    
    # Calculate terms of the Right-Hand Side (RHS)
    term_1 = 2**(1 - z)
    term_2 = f(z / 2)
    term_3 = f((z + 1) / 2)
    
    # Calculate RHS
    rhs = term_1 * term_2 * term_3
    
    print("The functional equation is: f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)")
    print("\nCalculated values for each term:")
    print(f"LHS: f(z)                = {lhs}")
    print(f"RHS part 1: 2^(1-z)       = {term_1}")
    print(f"RHS part 2: f(z/2)          = {term_2}")
    print(f"RHS part 3: f((z+1)/2)       = {term_3}")
    print(f"Total RHS                 = {rhs}")
    
    # Compare LHS and RHS
    if np.allclose(lhs, rhs):
        print("\nConclusion: The equality holds. The derived function is correct.\n")
    else:
        print("\nConclusion: The equality does not hold.\n")

# --- Main execution ---
# Choose a complex number for verification
z_test = 2.5 + 1.5j
verify_functional_equation(z_test)

# Verify the initial condition f(1) = sqrt(pi)
print("--- Verifying initial condition f(1) = sqrt(pi) ---")
f_1 = f(1)
sqrt_pi = np.sqrt(np.pi)
print(f"Calculated f(1) = {f_1}")
print(f"Expected sqrt(pi) = {sqrt_pi}")
if np.allclose(f_1, sqrt_pi):
    print("\nConclusion: The initial condition is satisfied.")
else:
    print("\nConclusion: The initial condition is NOT satisfied.")
