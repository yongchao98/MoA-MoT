import sympy

def solve_and_display():
    """
    This function defines the physical quantities as symbolic variables
    and prints the final expressions for the potential and electric field
    outside the sphere.
    """
    # Define the symbolic variables
    r, R, theta, E0, sigma1, sigma2 = sympy.symbols('r R theta E_0 sigma_1 sigma_2')
    
    # Define the potential outside the sphere (r > R) based on the derivation
    factor = (sigma1 - sigma2) / (sigma1 + 2 * sigma2)
    phi_out = -E0 * (r - factor * R**3 / r**2) * sympy.cos(theta)
    
    # Define the electric field components outside the sphere (r > R)
    E_r_out = E0 * (1 + 2 * factor * R**3 / r**3) * sympy.cos(theta)
    E_theta_out = -E0 * (1 - factor * R**3 / r**3) * sympy.sin(theta)

    # Print the final equations for the region r > R
    print("The question asks for the potential and electric field outside the sphere (r > R).")
    print("Based on the derivation, the correct expressions are:")
    print("-" * 50)
    
    print("Electric Potential Phi(r, theta) for r > R:")
    sympy.pprint(phi_out, use_unicode=True)
    print("\nWhich is equivalent to the expression in option B:")
    print("Phi(r, theta) = -E_0 * (r - (sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**2)) * cos(theta)\n")


    print("-" * 50)
    print("Electric Field E(r, theta) for r > R:")
    print("E_r component:")
    sympy.pprint(E_r_out, use_unicode=True)
    print("\nE_theta component:")
    sympy.pprint(E_theta_out, use_unicode=True)

    print("\nCombining these components gives the vector field E(r, theta), which matches option B.")
    print("E(r,theta) = E_0*[1 + 2*(sigma_1-sigma_2)*R^3/((sigma_1+2*sigma_2)*r^3)]*cos(theta) r_hat - E_0*[1 - (sigma_1-sigma_2)*R^3/((sigma_1+2*sigma_2)*r^3)]*sin(theta) theta_hat")

solve_and_display()