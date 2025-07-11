import numpy as np
from scipy.special import gamma

def verify_solution():
    """
    This function provides the derived solution to the functional equation
    and numerically verifies its correctness.
    """

    # --- Step 1: Explain the derived solution ---
    print("--- Derivation Summary ---")
    print("The functional equation is: f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)")
    print("The given condition is: f(1) = sqrt(pi)")
    print("By comparing with the Legendre duplication formula for the Gamma function,")
    print("we can deduce the solution's form.")
    print("\nThe explicit form of the function is f(z) = sqrt(pi) / Γ(z).\n")

    # --- Step 2: Define the derived function f(z) and print the final equation ---
    def f(z):
        return np.sqrt(np.pi) / gamma(z)

    print("--- Final Equation Details ---")
    print("f(z) = sqrt(pi) / Γ(z)")
    # Output the numbers in the final equation as requested
    print(f"The value of the constant pi is: {np.pi}")
    print(f"The value of the constant sqrt(pi) is: {np.sqrt(np.pi)}\n")

    # --- Step 3: Numerically verify the solution ---
    print("--- Numerical Verification ---")

    # Part A: Check the condition f(1) = sqrt(pi)
    z1 = 1
    f_at_1 = f(z1)
    sqrt_pi = np.sqrt(np.pi)
    print(f"Verifying the condition at z = {z1}:")
    print(f"Calculated f(1) = {f_at_1}")
    print(f"Expected sqrt(pi) = {sqrt_pi}")
    print(f"Is f(1) equal to sqrt(pi)? {np.isclose(f_at_1, sqrt_pi)}\n")

    # Part B: Check the functional equation for a sample complex number
    z_test = 2.5 + 1.5j
    print(f"Verifying the functional equation for a test value z = {z_test}:")

    # Calculate the Left-Hand Side (LHS)
    lhs = f(z_test)
    print(f"LHS = f(z) = {lhs}")

    # Calculate the Right-Hand Side (RHS)
    rhs = 2**(1 - z_test) * f(z_test / 2) * f((z_test + 1) / 2)
    print(f"RHS = 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs}")

    # Compare LHS and RHS
    print(f"Is LHS approximately equal to RHS? {np.isclose(lhs, rhs)}")

verify_solution()