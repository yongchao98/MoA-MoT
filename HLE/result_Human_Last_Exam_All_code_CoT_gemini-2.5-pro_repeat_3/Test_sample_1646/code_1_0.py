import cmath
import math
from scipy.special import gamma

def run_verification():
    """
    This script verifies that the function f(z) = sqrt(pi) / Gamma(z)
    is the solution to the given functional equation and conditions.
    """

    # The derived solution is f(z) = sqrt(pi) / Gamma(z).
    # We define this function in Python.
    def f(z):
        """
        Calculates sqrt(pi) / Gamma(z).
        """
        return cmath.sqrt(math.pi) / gamma(z)

    print("Step 1: Verifying the functional equation f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)")
    print("-" * 70)
    
    # We choose an arbitrary complex number to test the identity.
    z = 3.2 + 2.4j
    print(f"Testing for z = {z}\n")

    # Calculate the Left-Hand Side (LHS) of the equation
    lhs = f(z)
    
    # Calculate the Right-Hand Side (RHS) of the equation
    term1 = 2**(1 - z)
    term2 = f(z / 2)
    term3 = f((z + 1) / 2)
    rhs = term1 * term2 * term3

    print(f"LHS value: f(z) = {lhs}")
    print(f"RHS value: 2^(1 - z) * f(z/2) * f((z+1)/2) = {rhs}\n")
    
    if cmath.isclose(lhs, rhs):
        print("Result: LHS and RHS are approximately equal. The functional equation holds.\n")
    else:
        print("Result: LHS and RHS are not equal. The functional equation does not hold.\n")

    print("Step 2: Verifying the condition f(1) = sqrt(pi)")
    print("-" * 70)
    
    # Calculate f(1) using our function
    f_of_1 = f(1)
    
    # The expected value is sqrt(pi)
    sqrt_pi = math.sqrt(math.pi)
    
    print(f"Calculated value: f(1) = {f_of_1.real}")
    print(f"Expected value:   sqrt(pi) = {sqrt_pi}\n")
    
    if math.isclose(f_of_1.real, sqrt_pi):
        print("Result: The condition f(1) = sqrt(pi) is satisfied.\n")
    else:
        print("Result: The condition is not satisfied.\n")

    print("Conclusion:")
    print("The explicit form of the function is f(z) = c / Gamma(z), where c is a constant.")
    print(f"The constant c is sqrt(pi), with a value of approximately {sqrt_pi}.")
    print("So the final equation is: f(z) = sqrt(pi) / Gamma(z)")

run_verification()