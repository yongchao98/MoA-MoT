import numpy as np
from scipy.special import gamma

def solve_and_verify():
    """
    This function demonstrates the solution to the functional equation and verifies it numerically.
    
    The functional equation is: f(z) = 2**(1 - z) * f(z/2) * f((z+1)/2)
    The condition is: f(1) = sqrt(pi)

    The derivation shows that the solution is f(z) = sqrt(pi) / Gamma(z).
    We will verify this solution.
    """
    
    # Define the value of pi
    pi_val = np.pi

    # Define the proposed solution function f(z)
    def f(z):
        """Calculates sqrt(pi) / Gamma(z)"""
        return np.sqrt(pi_val) / gamma(z)

    print("Step 1: Verify the functional equation f(z) = 2**(1 - z) * f(z/2) * f((z+1)/2)")
    
    # We choose an arbitrary complex number for z to test the equation
    z_test = 3.5 + 2.5j

    # Calculate the left-hand side (LHS) of the equation
    lhs = f(z_test)

    # Calculate the components of the right-hand side (RHS)
    term1 = 2**(1 - z_test)
    term2 = f(z_test / 2)
    term3 = f((z_test + 1) / 2)
    rhs = term1 * term2 * term3
    
    print(f"Testing for z = {z_test}")
    print(f"LHS = f(z) = {lhs}")
    print(f"RHS = 2**(1-z) * f(z/2) * f((z+1)/2) = {rhs}")

    # Check if LHS and RHS are close enough to be considered equal
    if np.isclose(lhs, rhs):
        print("Verification successful: The functional equation holds.\n")
    else:
        print("Verification failed: The functional equation does not hold.\n")

    print("Step 2: Verify the condition f(1) = sqrt(pi)")

    # Calculate f(1) using the derived function
    f_of_1 = f(1)
    
    # The expected value is sqrt(pi)
    expected_value = np.sqrt(pi_val)
    
    print(f"Calculated value of f(1) is: {f_of_1}")
    print(f"Expected value of sqrt(pi) is: {expected_value}")
    
    if np.isclose(f_of_1, expected_value):
        print("Verification successful: The condition f(1) = sqrt(pi) holds.\n")
    else:
        print("Verification failed: The condition does not hold.\n")
        
    print("Conclusion: The explicit form of the function satisfying the given properties is:")
    print("f(z) = sqrt(pi) / Gamma(z)")
    
    # To satisfy the instruction "output each number in the final equation!",
    # we print the values of the constants involved in the final explicit form.
    print("\nConstants in the final equation f(z) = sqrt(pi) / Gamma(z):")
    print(f"pi = {pi_val}")
    print(f"sqrt(pi) = {np.sqrt(pi_val)}")


solve_and_verify()