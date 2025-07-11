import numpy as np
from scipy.special import gamma

def f(z):
    """
    Calculates the value of the derived function f(z) = sqrt(pi) / Gamma(z).
    """
    return np.sqrt(np.pi) / gamma(z)

def verify_solution():
    """
    Verifies the derived function against the given conditions.
    """
    # 1. Verify the initial condition f(1) = sqrt(pi)
    val_at_1 = f(1)
    expected_val_at_1 = np.sqrt(np.pi)
    print(f"Verification of f(1) = sqrt(pi):")
    print(f"  f(1) = {val_at_1:.6f}")
    print(f"  sqrt(pi) = {expected_val_at_1:.6f}")
    assert np.isclose(val_at_1, expected_val_at_1), "Verification failed for f(1)"
    print("  Condition f(1) = sqrt(pi) is satisfied.\n")

    # 2. Verify the functional equation for a sample complex number
    z = 2.5 + 1.5j
    
    # Calculate LHS of the functional equation: f(z)
    lhs = f(z)

    # Calculate RHS of the functional equation: 2^(1-z) * f(z/2) * f((z+1)/2)
    rhs = (2**(1 - z)) * f(z / 2) * f((z + 1) / 2)

    print(f"Verification of the functional equation for z = {z}:")
    print(f"  LHS = f(z) = {lhs:.6f}")
    print(f"  RHS = 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs:.6f}")
    assert np.isclose(lhs, rhs), "Functional equation verification failed."
    print("  Functional equation is satisfied.\n")
    
    # Print the final derived explicit form of the function
    pi_val = np.pi
    print("Conclusion:")
    print("The explicit form of the function is f(z) = sqrt(pi) / Gamma(z)")
    print(f"where pi is the constant {pi_val}")
    # The final equation involves the number 1 from the exponent (1-z) and the constant pi.
    print("Numbers in the final equation:")
    print(f"  - pi: {pi_val}")
    print(f"  - 1: 1")


if __name__ == "__main__":
    verify_solution()