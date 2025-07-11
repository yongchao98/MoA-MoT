import sympy as sp

def display_solution():
    """
    This function displays the derived expressions for the electric potential and
    electric field outside the sphere.
    """
    # Define the symbols used in the equations
    r, theta, R, E0, sigma1, sigma2 = sp.symbols('r, theta, R, E_0, sigma_1, sigma_2', real=True)

    # Expression for the potential outside the sphere (r > R)
    # This corresponds to the result in choice B
    potential_out = -E0 * (r - (sigma1 - sigma2) * R**3 / ((sigma1 + 2 * sigma2) * r**2)) * sp.cos(theta)
    
    # Calculate the electric field components from the potential using E = -grad(Phi)
    E_r_out = -sp.diff(potential_out, r)
    E_theta_out = -sp.diff(potential_out, theta) / r

    # Simplify the expressions to match the format in the answer choices
    E_r_out_simplified = E0 * (1 + 2 * (sigma1 - sigma2) * R**3 / ((sigma1 + 2 * sigma2) * r**3)) * sp.cos(theta)
    E_theta_out_simplified = -E0 * (1 - (sigma1 - sigma2) * R**3 / ((sigma1 + 2 * sigma2) * r**3)) * sp.sin(theta)

    print("The electric potential and electric field in the region outside the sphere (r > R) are:")
    
    print("\nPotential Phi(r, theta):")
    # To match the output format of the problem, we print a string representation
    phi_str = f"-E_0 * (r - (({sigma1} - {sigma2}) * R**3) / (({sigma1} + 2*{sigma2})*r**2)) * cos(theta)"
    print(phi_str)

    print("\nElectric Field E(r, theta):")
    Er_str = f"Er = E_0 * (1 + (2*({sigma1} - {sigma2}) * R**3) / (({sigma1} + 2*{sigma2})*r**3)) * cos(theta)"
    Etheta_str = f"E_theta = -E_0 * (1 - (({sigma1} - {sigma2}) * R**3) / (({sigma1} + 2*{sigma2})*r**3)) * sin(theta)"
    print(Er_str)
    print(Etheta_str)

    print("\nThese expressions correspond to the solution presented in option B.")

# Execute the function to print the solution
display_solution()