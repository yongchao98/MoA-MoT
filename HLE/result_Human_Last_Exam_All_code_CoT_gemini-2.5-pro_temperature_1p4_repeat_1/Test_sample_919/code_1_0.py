import sympy

def solve_emi_shielding_force():
    """
    This function uses the sympy library to symbolically construct the formula for the
    force per unit area on the conducting plane in the described EMI shielding application.
    It prints the components of the formula and the final result.
    """
    # Define the symbolic variables used in the problem.
    # We assume they are real and positive for simplicity.
    K0, a, y, d, mu, mu0 = sympy.symbols('K_0 a y d mu mu_0', real=True, positive=True)

    print("Deriving the force per unit area on the conductor at x = d.")
    print("The final formula has the form: F/area = (Numerator / Denominator) * i_x\n")

    # The numerator is derived from the magnetic pressure formula, proportional to B^2.
    # Numerator = mu_0 * (K_0 * sin(a*y))^2
    numerator_term = mu0 * K0**2 * sympy.sin(a*y)**2
    
    print("Component 1: Numerator of the force magnitude expression")
    print("="*50)
    sympy.pprint(numerator_term, use_unicode=True)
    print("="*50 + "\n")

    # The denominator arises from solving the boundary conditions for the magnetic fields.
    # The base of the denominator is the term [cosh(ad) + (mu_0/mu) * sinh(ad)]
    # The full denominator is 2 * [base]^2
    denominator_base = sympy.cosh(a*d) + (mu0 / mu) * sympy.sinh(a*d)
    denominator_term = 2 * denominator_base**2
    
    print("Component 2: Denominator of the force magnitude expression")
    print("="*50)
    sympy.pprint(denominator_term, use_unicode=True)
    print("="*50 + "\n")
    
    # The complete expression for the force magnitude
    force_magnitude = numerator_term / denominator_term
    
    # The force vector is in the x-direction (i_x)
    print("Final Result: Force per unit area vector")
    print("="*50)
    print("f/area =")
    sympy.pprint(force_magnitude, use_unicode=True)
    print("* i_x")
    print("="*50 + "\n")
    
    print("This symbolic expression corresponds to answer choice C.")

if __name__ == '__main__':
    solve_emi_shielding_force()