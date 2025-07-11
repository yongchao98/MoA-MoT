def compute_poynting_vector():
    """
    This function derives and prints the Poynting vector for the given scenario.
    It provides the symbolic expression for the vector both inside and outside the rod.

    Variables used in the expressions:
    R: radius of the cylindrical rod
    rho: uniform volume charge density of the rod
    E: magnitude of the uniform external electric field (along the z-axis)
    v: speed of the rod (along the z-axis)
    r: radial distance from the center of the rod
    epsilon_0: permittivity of free space
    mu_0: permeability of free space
    r_hat, z_hat: unit vectors in the radial and axial directions
    """

    print("This script computes the Poynting vector S = (1/mu_0) * (E_total x B).")
    print("The final result is presented in two parts: one for inside the rod (r < R) and one for outside the rod (r > R).")
    print("-" * 50)

    # --- Part 1: Inside the rod (r < R) ---
    print("For a point inside the rod (r < R):\n")

    # Electric Field components
    E_r_in = "(rho * r) / (2 * epsilon_0)"
    E_z_in = "E"
    
    # Magnetic Field component
    B_phi_in = "(mu_0 * rho * v * r) / 2"

    print("The total electric field has components:")
    print(f"  E_r = {E_r_in}")
    print(f"  E_z = {E_z_in}\n")

    print("The magnetic field has component:")
    print(f"  B_phi = {B_phi_in}\n")

    # Poynting Vector components
    # S = (1/mu_0) * (E_r * B_phi * z_hat - E_z * B_phi * r_hat)
    S_r_in_final = "-(E * rho * v * r) / 2"
    S_z_in_final = "(rho**2 * v * r**2) / (4 * epsilon_0)"

    print("The resulting Poynting vector S = S_r * r_hat + S_z * z_hat has components:")
    print("  S_r (radial component) = " + S_r_in_final)
    print("  S_z (axial component)  = " + S_z_in_final)
    print("-" * 50)


    # --- Part 2: Outside the rod (r > R) ---
    print("For a point outside the rod (r > R):\n")

    # Electric Field components
    E_r_out = "(rho * R**2) / (2 * epsilon_0 * r)"
    E_z_out = "E"
    
    # Magnetic Field component
    B_phi_out = "(mu_0 * rho * v * R**2) / (2 * r)"

    print("The total electric field has components:")
    print(f"  E_r = {E_r_out}")
    print(f"  E_z = {E_z_out}\n")

    print("The magnetic field has component:")
    print(f"  B_phi = {B_phi_out}\n")
    
    # Poynting Vector components
    S_r_out_final = "-(E * rho * v * R**2) / (2 * r)"
    S_z_out_final = "(rho**2 * v * R**4) / (4 * epsilon_0 * r**2)"
    
    print("The resulting Poynting vector S = S_r * r_hat + S_z * z_hat has components:")
    print("  S_r (radial component) = " + S_r_out_final)
    print("  S_z (axial component)  = " + S_z_out_final)
    print("-" * 50)
    
    print("\nFinal symbolic equations for the Poynting vector S:")
    print("Inside the rod (r < R):")
    final_eq_in = f"S = ({S_r_in_final}) r_hat + ({S_z_in_final}) z_hat"
    print(final_eq_in)

    print("\nOutside the rod (r > R):")
    final_eq_out = f"S = ({S_r_out_final}) r_hat + ({S_z_out_final}) z_hat"
    print(final_eq_out)


if __name__ == "__main__":
    compute_poynting_vector()