import sympy
from sympy import init_printing

def solve_electrodynamics_problem():
    """
    This function defines and prints the expressions for the electric potential and
    electric field outside a conducting sphere in a conducting medium, under an
    external uniform electric field in a steady state.
    """
    # Initialize pretty printing for mathematical formulas
    init_printing(use_unicode=True)

    # Define the symbolic variables
    r, R, theta = sympy.symbols('r R theta')
    E0, sigma1, sigma2 = sympy.symbols('E_0 sigma_1 sigma_2')
    
    # Define a shorthand for the coefficient term
    C = (sigma1 - sigma2) / (sigma1 + 2 * sigma2)

    # Define the electric potential Phi for r > R (outside the sphere)
    # Phi(r, theta) = -E0 * (r - C * R**3 / r**2) * cos(theta)
    phi_out = -E0 * (r - ( (sigma1 - sigma2) * R**3 ) / ( (sigma1 + 2 * sigma2) * r**2 ) ) * sympy.cos(theta)
    
    # Calculate the components of the electric field E = -grad(Phi)
    # E_r = -d(Phi)/dr
    E_r_out = -sympy.diff(phi_out, r)
    
    # E_theta = -1/r * d(Phi)/d(theta)
    E_theta_out = - (1/r) * sympy.diff(phi_out, theta)

    # Print the results for the region outside the sphere (r > R)
    print("The electric potential Phi(r, theta) outside the sphere (r > R) is:")
    sympy.pprint(phi_out)
    
    print("\n" + "="*50 + "\n")
    
    print("The electric field E(r, theta) outside the sphere (r > R) has components:")
    print("\nRadial component (E_r):")
    # Simplify and display E_r
    sympy.pprint(sympy.simplify(E_r_out))
    
    print("\nAngular component (E_theta):")
    # Simplify and display E_theta
    sympy.pprint(sympy.simplify(E_theta_out))

    print("\n" + "="*50 + "\n")
    print("The full vector for the electric field E(r, theta) outside the sphere is E_r * r_hat + E_theta * theta_hat.")
    print("Based on the derivation, the correct option is B.")

solve_electrodynamics_problem()