import cmath
import math

def f(z):
    """
    Calculates the derived function f(z) = sqrt(pi) / Gamma(z).
    """
    # Using cmath.gamma for complex number support
    return math.sqrt(math.pi) / cmath.gamma(z)

def verify_functional_equation(z):
    """
    Verifies the functional equation for a given complex number z.
    The equation is f(z) = 2^(1-z) * f(z/2) * f((z+1)/2).
    """
    print(f"Verifying for z = {z}\n")

    # --- Calculate LHS ---
    lhs = f(z)

    # --- Calculate RHS components ---
    term1 = 2**(1 - z)
    term2 = f(z / 2)
    term3 = f((z + 1) / 2)

    # --- Calculate full RHS ---
    rhs = term1 * term2 * term3

    # --- Print the results ---
    print("The functional equation is: f(z) = 2**(1-z) * f(z/2) * f((z+1)/2)\n")
    print("--- Calculated Values ---")
    print(f"LHS: f(z) = {lhs}")
    print(f"RHS: {rhs}\n")

    print("--- Breakdown of the RHS Calculation ---")
    # Using the requirement to output each number in the final equation
    print(f"  2**(1-z)       = {term1}")
    print(f"  f(z/2)         = {term2}")
    print(f"  f((z+1)/2)     = {term3}")
    print(f"  Product (RHS)  = {term1 * term2 * term3}")
    
    # --- Check if they are close enough ---
    print("\n--- Verification ---")
    if cmath.isclose(lhs, rhs):
        print("LHS and RHS are equal. The solution is correct.")
    else:
        print("LHS and RHS are not equal. There might be an issue.")


# Choose a complex number for verification
z_to_test = 2.5 + 1.5j
verify_functional_equation(z_to_test)