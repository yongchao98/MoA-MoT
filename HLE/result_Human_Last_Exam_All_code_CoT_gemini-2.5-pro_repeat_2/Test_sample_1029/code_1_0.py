import sympy

def compute_poynting_vector():
    """
    This function derives and prints the symbolic expression for the Poynting vector
    for a moving charged rod in an external electric field.

    The final result is presented for two regions: inside the rod (r < R)
    and outside the rod (r > R).

    Nomenclature:
    R: radius of the rod
    rho: volume charge density
    E: magnitude of the external electric field (along the z-axis)
    v: speed of the rod (along the z-axis)
    r: radial distance from the center of the rod
    epsilon_0: permittivity of free space
    mu_0: permeability of free space
    r_hat, phi_hat, k_hat: unit vectors in cylindrical coordinates
    """

    # --- Derivations (for clarity, not used in final print) ---
    # E_ext = E * k_hat
    # E_rod_in = (rho * r / (2 * epsilon_0)) * r_hat
    # E_total_in = E_rod_in + E_ext
    #
    # B_in = (mu_0 * rho * v * r / 2) * phi_hat
    #
    # S_in = (1/mu_0) * cross(E_total_in, B_in)
    #
    # E_rod_out = (rho * R**2 / (2 * epsilon_0 * r)) * r_hat
    # E_total_out = E_rod_out + E_ext
    #
    # B_out = (mu_0 * rho * v * R**2 / (2*r)) * phi_hat
    #
    # S_out = (1/mu_0) * cross(E_total_out, B_out)
    # --- End Derivations ---

    print("The Poynting vector S = (1/mu_0) * (E_total x B_total) is calculated for the system.")
    print("The result is expressed in cylindrical coordinates (r, phi, z), with unit vectors r_hat, phi_hat, and k_hat (along the axis).")
    print("-" * 70)

    # --- Inside the rod (r < R) ---
    print("For the region INSIDE the rod (r < R):")
    
    # The components of the Poynting vector
    Sr_in_coeff_1 = "- (E * rho * v)"
    Sr_in_coeff_2 = "2"
    Sr_in_var = "r"
    
    Sk_in_coeff_1 = "(rho**2 * v)"
    Sk_in_coeff_2 = "(4 * epsilon_0)"
    Sk_in_var = "r**2"

    print("The radial component of the Poynting vector (S_r * r_hat) is:")
    print(f"  Coefficient: {Sr_in_coeff_1} / {Sr_in_coeff_2}")
    print(f"  Variable part: {Sr_in_var}")
    print(f"  S_r = ({Sr_in_coeff_1} * {Sr_in_var}) / {Sr_in_coeff_2}")
    
    print("\nThe axial component of the Poynting vector (S_k * k_hat) is:")
    print(f"  Coefficient: {Sk_in_coeff_1} / {Sk_in_coeff_2}")
    print(f"  Variable part: {Sk_in_var}")
    print(f"  S_k = ({Sk_in_coeff_1} * {Sk_in_var}) / {Sk_in_coeff_2}")

    print("\nTherefore, the full Poynting vector inside the rod is:")
    print(f"  S_inside = - (E * rho * v * r / 2) r_hat + (rho**2 * v * r**2 / (4 * epsilon_0)) k_hat")
    print("-" * 70)
    
    # --- Outside the rod (r > R) ---
    print("For the region OUTSIDE the rod (r > R):")

    Sr_out_coeff_1 = "- (E * rho * v * R**2)"
    Sr_out_coeff_2 = "2"
    Sr_out_var = "1/r"

    Sk_out_coeff_1 = "(rho**2 * v * R**4)"
    Sk_out_coeff_2 = "(4 * epsilon_0)"
    Sk_out_var = "1/r**2"

    print("The radial component of the Poynting vector (S_r * r_hat) is:")
    print(f"  Coefficient: {Sr_out_coeff_1} / {Sr_out_coeff_2}")
    print(f"  Variable part: {Sr_out_var}")
    print(f"  S_r = ({Sr_out_coeff_1}) / ({Sr_out_coeff_2} * r)")

    print("\nThe axial component of the Poynting vector (S_k * k_hat) is:")
    print(f"  Coefficient: {Sk_out_coeff_1} / {Sk_out_coeff_2}")
    print(f"  Variable part: {Sk_out_var}")
    print(f"  S_k = ({Sk_out_coeff_1}) / ({Sk_out_coeff_2} * r**2)")
    
    print("\nTherefore, the full Poynting vector outside the rod is:")
    print(f"  S_outside = - (E * rho * v * R**2 / (2 * r)) r_hat + (rho**2 * v * R**4 / (4 * epsilon_0 * r**2)) k_hat")

if __name__ == '__main__':
    compute_poynting_vector()