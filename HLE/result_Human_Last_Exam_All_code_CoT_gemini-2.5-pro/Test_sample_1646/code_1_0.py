import cmath
import math
from scipy.special import gamma

def solve_and_verify():
    """
    This function states the solution to the functional equation and verifies it numerically.
    """
    
    solution_formula = "f(z) = sqrt(pi) / Gamma(z)"
    print(f"The explicit form of f(z) is determined to be: {solution_formula}\n")
    
    # Define the function based on the derived solution
    def f(z):
        # The Gamma function has poles at z = 0, -1, -2, ...
        # At these points, its reciprocal 1/Gamma(z) is 0. Since f(z) must be
        # analytic everywhere, this is the correct behavior.
        if z.imag == 0 and z.real <= 0 and z.real == int(z.real):
            return 0.0
        return math.sqrt(math.pi) / gamma(z)

    print("--- Verification of the Solution ---")

    # 1. Verify the initial condition f(1) = sqrt(pi)
    z1 = 1
    f_at_1 = f(complex(z1))
    sqrt_pi = math.sqrt(math.pi)
    
    print(f"\n1. Verifying the condition f(1) = sqrt(pi):")
    print(f"   The constant pi is approximately {math.pi}")
    print(f"   Calculated f(1) = {f_at_1.real}")
    print(f"   Expected value sqrt(pi) = {sqrt_pi}")
    if math.isclose(f_at_1.real, sqrt_pi):
        print("   Condition is satisfied.")
    else:
        print("   Condition is NOT satisfied.")

    # 2. Verify the functional equation for a sample complex value
    # The equation is: f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)
    z_test = complex(2.5, 1.5)
    
    lhs = f(z_test)
    
    c_factor = 2**(1 - z_test)
    term1 = f(z_test / 2)
    term2 = f((z_test + 1) / 2)
    rhs = c_factor * term1 * term2
    
    print(f"\n2. Verifying the functional equation for z = {z_test}:")
    print(f"   The numbers in the equation are: 2, 1")
    print(f"   LHS: f({z_test}) = {lhs}")
    print(f"   RHS: 2^(1 - ({z_test})) * f({z_test / 2}) * f({(z_test + 1) / 2})")
    print(f"        = {c_factor} * {term1} * {term2}")
    print(f"        = {rhs}")

    if cmath.isclose(lhs, rhs, rel_tol=1e-9):
        print("   Functional equation is satisfied.")
    else:
        print("   Functional equation is NOT satisfied.")

# Execute the verification
solve_and_verify()