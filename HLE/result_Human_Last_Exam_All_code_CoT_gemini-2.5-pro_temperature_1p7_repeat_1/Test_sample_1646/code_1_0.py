import numpy as np
from scipy.special import gamma

def solve_functional_equation():
    """
    This script determines and verifies the explicit form of the function f(z).
    
    The function f(z) is analytic and satisfies:
    1. f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)
    2. f(1) = sqrt(pi)

    Through mathematical derivation using the Legendre duplication formula for the
    Gamma function, f(z) is determined to be sqrt(pi) / Gamma(z).
    """
    
    # 1. State the final derived form and the constants involved.
    print("The explicit form of the function is f(z) = sqrt(pi) / Gamma(z).")
    
    pi_val = np.pi
    sqrt_pi_val = np.sqrt(pi_val)
    
    print("\nConstants in the final equation:")
    print(f"The value of pi is: {pi_val}")
    print(f"The value of sqrt(pi) is: {sqrt_pi_val}")
    
    # 2. Define the derived function f(z)
    def f(z):
        return np.sqrt(np.pi) / gamma(z)

    # 3. Numerically verify the functional equation for a sample complex value
    print("\n--- Verification ---")
    z_test = 2.5 + 1.5j
    
    lhs = f(z_test)
    rhs = (2**(1 - z_test)) * f(z_test / 2) * f((z_test + 1) / 2)
    
    print(f"Testing the functional equation with z = {z_test}")
    print(f"LHS f(z)                   = {lhs}")
    print(f"RHS 2^(1-z)f(z/2)f((z+1)/2) = {rhs}")
    
    # Using np.isclose for robust floating point comparison
    if np.isclose(lhs, rhs):
        print("Functional equation holds true.")
    else:
        print("Functional equation does not hold.")

    # 4. Numerically verify the given condition f(1) = sqrt(pi)
    f_at_1 = f(1)
    
    print("\nTesting the initial condition f(1) = sqrt(pi)")
    print(f"Calculated f(1) = {f_at_1}")
    print(f"Expected sqrt(pi) = {sqrt_pi_val}")
    
    if np.isclose(f_at_1, sqrt_pi_val):
        print("Initial condition holds true.")
    else:
        print("Initial condition does not hold.")
        
    print("\n--- Conclusion ---")
    print("The explicit form of f(z) is sqrt(pi) / Gamma(z).")

solve_functional_equation()
