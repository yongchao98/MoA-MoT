import math

def solve_and_verify():
    """
    This script verifies the derived solution for the functional equation and
    then prints the final explicit formula for f(z).
    """

    # The problem asks for the explicit form of f(z).
    # Based on the derivation, the solution is f(z) = sqrt(pi) / Gamma(z).

    # This script will verify the solution for a sample real number z > 0
    # using the standard `math` library.
    
    print("--- Verification of the Solution f(z) = sqrt(pi) / Gamma(z) ---")
    
    # Let's pick a sample real value for z
    z = 2.5
    
    # Define the function f(z) based on the derived formula
    f = lambda val: math.sqrt(math.pi) / math.gamma(val)
    
    # --- Verify the functional equation ---
    # Left-hand side: f(z)
    lhs = f(z)
    
    # Right-hand side: 2**(1-z) * f(z/2) * f((z+1)/2)
    rhs = (2**(1 - z)) * f(z / 2) * f((z + 1) / 2)
    
    print(f"\nVerifying the functional equation for z = {z}:")
    print(f"LHS = f(z) = {lhs}")
    print(f"RHS = 2**(1-z) * f(z/2) * f((z+1)/2) = {rhs}")
    
    if math.isclose(lhs, rhs):
        print("Result: LHS and RHS are approximately equal. The functional equation holds.")
    else:
        print("Result: LHS and RHS are not equal. Verification failed.")
        
    # --- Verify the initial condition ---
    f1 = f(1)
    sqrt_pi = math.sqrt(math.pi)
    
    print(f"\nVerifying the initial condition f(1) = sqrt(pi):")
    print(f"Calculated f(1) = {f1}")
    print(f"Expected value sqrt(pi) = {sqrt_pi}")
    
    if math.isclose(f1, sqrt_pi):
        print("Result: The initial condition is satisfied.")
    else:
        print("Result: The initial condition is not satisfied. Verification failed.")

    # --- Print the final explicit formula ---
    # The equation contains the number pi.
    print("\n--- Final Answer ---")
    print("The explicit form of the function f(z) is:")
    
    final_equation = "f(z) = sqrt(pi) / Gamma(z)"
    print(final_equation)

solve_and_verify()