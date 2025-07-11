import sympy as sp

def compute_poynting_vector():
    """
    This script calculates the Poynting vector S for a moving, charged cylindrical rod
    in a uniform external electric field.
    """
    
    # Define the symbolic variables used in the problem statement.
    # Unicode characters are used for better readability in the output.
    r, R, rho, E, v, epsilon_0, mu_0 = sp.symbols('r R ρ E v ε₀ μ₀', real=True, positive=True)

    print("This script computes the Poynting vector S for the given scenario.")
    print("The vector S has a radial component (S_r) and an axial component (S_z).")
    print("The final expressions depend on the radial distance 'r' from the axis.\n")
    
    # --- STEP 1: Define Fields Inside the Rod (r < R) ---

    # Electric Field (E = E_r * r_hat + E_z * z_hat)
    # E_r is from the rod's own charge (Gauss's Law)
    E_r_in = (rho * r) / (2 * epsilon_0)
    # E_z is the external uniform field
    E_z_in = E

    # Magnetic Field (B = B_phi * phi_hat)
    # B_phi is from the current created by the moving rod (Ampere's Law)
    B_phi_in = (mu_0 * rho * v * r) / 2

    # --- STEP 2: Calculate Poynting Vector Inside the Rod (r < R) ---

    # S = (1/mu_0) * (E x B)
    # The cross product E x B in cylindrical coordinates gives:
    # (E x B)_r = -E_z * B_phi
    # (E x B)_z = E_r * B_phi
    
    S_r_in = (1 / mu_0) * (-E_z_in * B_phi_in)
    S_z_in = (1 / mu_0) * (E_r_in * B_phi_in)

    # --- STEP 3: Define Fields Outside the Rod (r > R) ---

    # Electric Field
    E_r_out = (rho * R**2) / (2 * epsilon_0 * r)
    E_z_out = E
    
    # Magnetic Field
    B_phi_out = (mu_0 * rho * v * R**2) / (2 * r)

    # --- STEP 4: Calculate Poynting Vector Outside the Rod (r > R) ---
    
    S_r_out = (1 / mu_0) * (-E_z_out * B_phi_out)
    S_z_out = (1 / mu_0) * (E_r_out * B_phi_out)
    
    # --- STEP 5: Print the Final Equations ---
    
    # Using sp.pretty() for a more readable mathematical format
    print("--- For inside the rod (r < R): ---")
    print(f"The Poynting vector S(r) is given by: S = S_r * r_hat + S_z * z_hat")
    print("\nThe radial component S_r is:")
    sp.pretty_print(S_r_in)
    print("\nThe axial component S_z is:")
    sp.pretty_print(S_z_in)
    
    print("\n" + "="*40 + "\n")
    
    print("--- For outside the rod (r > R): ---")
    print(f"The Poynting vector S(r) is given by: S = S_r * r_hat + S_z * z_hat")
    print("\nThe radial component S_r is:")
    sp.pretty_print(S_r_out)
    print("\nThe axial component S_z is:")
    sp.pretty_print(S_z_out)


compute_poynting_vector()