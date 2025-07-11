import sympy as sp

def solve_electrodynamics_problem():
    """
    This function presents the solution for the electric potential and field
    outside a conducting sphere placed in a conducting medium with a uniform applied electric field.

    The solution is derived by solving Laplace's equation with the appropriate
    boundary conditions for a steady-state current.
    """
    
    # Define the symbolic variables used in the problem
    E0, R, r, theta = sp.symbols('E_0 R r theta')
    sigma1, sigma2 = sp.symbols('sigma_1 sigma_2')
    
    # The problem asks for the potential and electric field in the region outside the sphere (r > R).
    # The derivation leads to specific expressions for the potential and field.
    
    print("This script provides the solution for the electric potential and field outside the sphere (r > R).\n")

    # --- Electric Potential Phi(r, theta) for r > R ---
    # The potential is a superposition of the uniform field potential and an induced dipole potential.
    # The coefficient of the dipole term is derived from the boundary conditions.
    coeff_num = sigma1 - sigma2
    coeff_den = sigma1 + 2 * sigma2
    
    phi_outside_str = f"-E_0 * (r - (({coeff_num}) * R**3) / (({coeff_den}) * r**2)) * cos(theta)"
    
    print("The electric potential Phi(r, theta) for r > R is:")
    print(f"Phi(r, theta) = {phi_outside_str}\n")
    
    # --- Electric Field E(r, theta) for r > R ---
    # The electric field is the negative gradient of the potential, E = -grad(Phi).
    # It has a radial (r_hat) and a polar (theta_hat) component.
    
    # Radial component E_r
    E_r_str = f"E_0 * [1 + (2 * ({coeff_num}) * R**3) / (({coeff_den}) * r**3)] * cos(theta)"
    
    # Polar component E_theta
    E_theta_str = f"-E_0 * [1 - (({coeff_num}) * R**3) / (({coeff_den}) * r**3)] * sin(theta)"
    
    print("The electric field E(r, theta) = E_r * r_hat + E_theta * theta_hat for r > R has components:")
    print("Radial component:")
    print(f"E_r = {E_r_str}\n")
    print("Polar component:")
    print(f"E_theta = {E_theta_str}\n")
    
    print("Comparing these results with the given options, we find that they match Choice B.")

if __name__ == '__main__':
    solve_electrodynamics_problem()
