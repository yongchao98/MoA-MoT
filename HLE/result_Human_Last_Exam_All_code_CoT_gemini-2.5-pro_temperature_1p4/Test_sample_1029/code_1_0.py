import sympy as sp

def compute_poynting_vector():
    """
    This function derives and prints the Poynting vector for the given setup.
    It does not take numerical inputs but rather shows the symbolic derivation.
    """
    # Define symbolic variables for clarity in the derivation steps
    E, rho, v, r, R, mu_0, epsilon_0 = sp.symbols('E rho v r R mu_0 epsilon_0', real=True, positive=True)

    # --- Step 1: Electric Field ---
    # External field
    E_ext_z = E
    # Field from the rod (radial component)
    E_rod_r_inside = (rho * r) / (2 * epsilon_0)
    E_rod_r_outside = (rho * R**2) / (2 * epsilon_0 * r)

    # Total electric field components
    E_total_r_inside = E_rod_r_inside
    E_total_z_inside = E_ext_z

    E_total_r_outside = E_rod_r_outside
    E_total_z_outside = E_ext_z

    # --- Step 2: Magnetic Field ---
    # Magnetic field is purely azimuthal (phi-direction)
    B_phi_inside = (mu_0 * rho * v * r) / 2
    B_phi_outside = (mu_0 * rho * v * R**2) / (2 * r)

    # --- Step 3: Cross Product E x B ---
    # In cylindrical coordinates, (E_r*r_hat + E_z*z_hat) x (B_phi*phi_hat)
    # = E_r * B_phi * (r_hat x phi_hat) + E_z * B_phi * (z_hat x phi_hat)
    # = E_r * B_phi * z_hat - E_z * B_phi * r_hat

    # Cross product components inside the rod (r <= R)
    ExB_r_inside = -E_total_z_inside * B_phi_inside
    ExB_z_inside = E_total_r_inside * B_phi_inside

    # Cross product components outside the rod (r > R)
    ExB_r_outside = -E_total_z_outside * B_phi_outside
    ExB_z_outside = E_total_r_outside * B_phi_outside

    # --- Step 4: Poynting Vector S = (1/mu_0) * (E x B) ---
    S_r_inside = ExB_r_inside / mu_0
    S_z_inside = ExB_z_inside / mu_0

    S_r_outside = ExB_r_outside / mu_0
    S_z_outside = ExB_z_outside / mu_0
    
    # --- Step 5: Print the results ---
    print("The Poynting vector S is calculated in two regions:")
    print("-" * 50)
    
    print("Inside the rod (r <= R):")
    # Using sp.pretty to format the output
    print("S_r (radial component) = ")
    sp.pprint(S_r_inside, use_unicode=False)
    print("\nS_z (axial component) = ")
    sp.pprint(S_z_inside, use_unicode=False)
    print(f"\nSo, S_inside = ({S_r_inside}) r_hat + ({S_z_inside}) z_hat")
    
    print("\n" + "-" * 50)
    
    print("Outside the rod (r > R):")
    print("S_r (radial component) = ")
    sp.pprint(S_r_outside, use_unicode=False)
    print("\nS_z (axial component) = ")
    sp.pprint(S_z_outside, use_unicode=False)
    print(f"\nSo, S_outside = ({S_r_outside}) r_hat + ({S_z_outside}) z_hat")


if __name__ == '__main__':
    compute_poynting_vector()