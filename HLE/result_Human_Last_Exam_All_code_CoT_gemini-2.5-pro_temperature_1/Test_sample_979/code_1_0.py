import sympy as sp

def display_magnetic_field_solution():
    """
    This function displays the symbolic solution for the magnetic field H
    inside and outside a spherical shell with a given surface current.
    """

    # Define the symbolic variables
    mu, mu0, K0, R, r, theta = sp.symbols('mu mu_0 K_0 R r theta', real=True, positive=True)
    
    # Define unit vectors symbolically (as strings for display)
    r_hat = "r_hat"
    theta_hat = "theta_hat"
    z_hat = "z_hat"

    # Common denominator term from the derivation
    denominator = 1 + (2 * mu0) / mu

    # --- Field Inside the Sphere (r < R) ---
    # Coefficient for the inside field
    h_in_coeff = (2 * mu0 / mu) * (K0 / denominator)
    
    # The field inside is uniform and in the z-direction
    # H_in = h_in_coeff * z_hat
    
    print("Derived Magnetic Field Solution:\n")
    print("Region inside the sphere (r < R):")
    # We construct the string representation to match the format in the answer choices
    h_in_str = f"H_in = ( (2 * mu0 / mu) * K0 / (1 + 2 * mu0 / mu) ) * {z_hat}"
    print(h_in_str)
    
    # You can also see the simplified form:
    # print(f"Simplified H_in = {sp.simplify(h_in_coeff)} * {z_hat}")
    
    print("\n" + "="*40 + "\n")
    
    # --- Field Outside the Sphere (r > R) ---
    # Coefficient for the outside field
    h_out_coeff = K0 / denominator
    
    # The field outside is a dipole field
    # H_out = h_out_coeff * (R**3 / r**3) * (2*sp.cos(theta)*r_hat + sp.sin(theta)*theta_hat)
    
    print("Region outside the sphere (r > R):")
    # Constructing the string representation for the outside field
    h_out_str = (f"H_out = ( K0 / (1 + 2 * mu0 / mu) ) * (R**3 / r**3) * "
                 f"(2*cos({theta})*{r_hat} + sin({theta})*{theta_hat})")
    print(h_out_str)
    
    # You can also see the simplified form:
    # print(f"Simplified H_out = {sp.simplify(h_out_coeff)} * (R**3/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")

if __name__ == '__main__':
    display_magnetic_field_solution()
