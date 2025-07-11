import sympy
from sympy import I, pi, gamma, expand_func, simplify, Eq, solve

def solve_functional_equation():
    """
    This script verifies the solution to the given functional equation.
    
    The problem is to find f(z) given:
    1. f(z) = 2**(1 - z) * f(z/2) * f((z+1)/2)
    2. f(1) = sqrt(pi)
    
    The proposed solution is f(z) = sqrt(pi) / Gamma(z).
    We will verify this using symbolic mathematics.
    """
    
    # Define the symbolic variable z
    z = sympy.Symbol('z')

    # Define our proposed solution for f(z)
    # f(z) = sqrt(pi) / Gamma(z)
    f_z = sympy.sqrt(pi) / gamma(z)

    print("Proposed solution: f(z) = ", f_z)
    print("-" * 30)

    # Step 1: Verify the condition f(1) = sqrt(pi)
    print("Step 1: Verifying the initial condition f(1) = sqrt(pi)")
    
    # Substitute z = 1 into our function
    f_at_1 = f_z.subs(z, 1)
    
    # The gamma function at 1 is 1
    gamma_at_1 = gamma(1)
    
    print(f"f(1) = sqrt(pi) / Gamma(1)")
    print(f"Since Gamma(1) = {gamma_at_1},")
    print(f"f(1) = sqrt(pi) / {gamma_at_1} = {f_at_1}")
    
    # Check if f(1) equals sqrt(pi)
    if f_at_1 == sympy.sqrt(pi):
        print("The condition f(1) = sqrt(pi) is satisfied.")
    else:
        print("The condition f(1) = sqrt(pi) is NOT satisfied.")
    print("-" * 30)

    # Step 2: Verify the functional equation
    # f(z) = 2**(1 - z) * f(z/2) * f((z+1)/2)
    print("Step 2: Verifying the functional equation.")

    # Left-Hand Side (LHS) of the equation is f(z)
    lhs = f_z
    print("LHS of the equation: f(z) =", lhs)

    # Right-Hand Side (RHS) of the equation is 2**(1-z) * f(z/2) * f((z+1)/2)
    # Substitute f(z) with our proposed solution
    f_z_div_2 = f_z.subs(z, z/2)
    f_z_plus_1_div_2 = f_z.subs(z, (z+1)/2)
    
    rhs_expression_str = f"2**(1 - z) * ({f_z_div_2}) * ({f_z_plus_1_div_2})"
    rhs = 2**(1 - z) * f_z_div_2 * f_z_plus_1_div_2
    
    print("RHS of the equation:", rhs)

    # Simplify the RHS using the Legendre duplication formula,
    # which sympy's expand_func knows for gamma functions.
    # The formula is: Gamma(z/2)*Gamma((z+1)/2) = 2**(1-z)*sqrt(pi)*Gamma(z)
    
    simplified_rhs = simplify(rhs)
    
    print("\nSimplifying the RHS:")
    print("The denominator contains Gamma(z/2) * Gamma((z+1)/2).")
    print("Using Legendre's duplication formula, this simplifies.")
    print("Simplified RHS:", simplified_rhs)
    
    # Final check
    if lhs == simplified_rhs:
        print("\nConclusion: The LHS and simplified RHS are equal.")
        print("The function f(z) = sqrt(pi) / Gamma(z) is the correct solution.")
    else:
        print("\nConclusion: The function does not satisfy the equation.")

    # Final explicit form as requested
    print("\n" + "="*40)
    print("The explicit form of the function is:")
    # "output each number in the final equation" can be interpreted as printing constants
    numerator = sympy.sqrt(pi)
    denominator_str = "Gamma(z)"
    print(f"f(z) = {numerator} / {denominator_str}")
    print("="*40)

if __name__ == '__main__':
    solve_functional_equation()
