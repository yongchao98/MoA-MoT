import sympy

def compute_poynting_vector():
    """
    This script symbolically computes the Poynting vector for a moving, charged cylindrical rod
    in an external electric field.

    The Poynting vector S is given by: S = (1/mu_0) * E_total x B
    
    The calculation is performed for two regions: inside the rod (r < R) and outside the rod (r > R).
    """
    
    # 1. Define symbolic variables for the physical quantities.
    # r: radial distance from the center of the rod
    # R: radius of the rod
    # rho: volume charge density
    # v: speed of the rod along its axis
    # E: magnitude of the external electric field
    # mu_0: permeability of free space
    # epsilon_0: permittivity of free space
    r, R, rho, v, E, mu_0, epsilon_0 = sympy.symbols('r R ρ v E μ₀ ε₀', real=True, positive=True)

    print("Computing the Poynting vector S = (1/μ₀) * E x B for a moving charged rod.\n")

    # --- Case 1: Inside the rod (r < R) ---

    print("--- Inside the rod (r < R) ---")

    # E_total = E_rod + E_ext
    # E_rod is radial. From Gauss's Law: E_rod * 2*pi*r*L = (rho * pi*r^2*L) / epsilon_0
    E_rod_in = (rho * r) / (2 * epsilon_0)
    # The external field E is axial (z-direction).
    # The total electric field has a radial component (E_r) and an axial component (E_z).
    E_r_in = E_rod_in
    E_z_in = E

    # B is azimuthal (phi-direction). From Ampere's Law: B * 2*pi*r = mu_0 * I_enc = mu_0 * (J * pi*r^2)
    # where current density J = rho * v.
    B_phi_in = (mu_0 * rho * v * r) / 2

    # Calculate the components of the Poynting vector S = (1/mu_0) * (E x B)
    # In cylindrical coordinates, E x B = (E_r*r_hat + E_z*z_hat) x (B_phi*phi_hat)
    # = E_r*B_phi * (r_hat x phi_hat) + E_z*B_phi * (z_hat x phi_hat)
    # = (E_r*B_phi * z_hat) - (E_z*B_phi * r_hat)
    # So, S_r = -E_z*B_phi and S_z = E_r*B_phi.
    
    S_r_in = -E_z_in * B_phi_in / mu_0
    S_z_in = E_r_in * B_phi_in / mu_0

    # Print results for inside the rod, showing each component of the vector.
    print("The Poynting vector S has a radial component (S_r) and an axial component (S_z).")
    print("\nRadial component S_r:")
    sympy.pprint(S_r_in, use_unicode=True)
    print("\nAxial component S_z:")
    sympy.pprint(S_z_in, use_unicode=True)
    print("-" * 40)


    # --- Case 2: Outside the rod (r > R) ---

    print("\n--- Outside the rod (r > R) ---")

    # E_rod from Gauss's Law: E_rod * 2*pi*r*L = (rho * pi*R^2*L) / epsilon_0
    E_rod_out = (rho * R**2) / (2 * epsilon_0 * r)
    E_r_out = E_rod_out
    E_z_out = E

    # B from Ampere's Law: B * 2*pi*r = mu_0 * I_enc = mu_0 * (rho*v * pi*R^2)
    B_phi_out = (mu_0 * rho * v * R**2) / (2 * r)

    # Calculate the components of the Poynting vector S.
    S_r_out = -E_z_out * B_phi_out / mu_0
    S_z_out = E_r_out * B_phi_out / mu_0

    # Print results for outside the rod.
    print("The Poynting vector S has a radial component (S_r) and an axial component (S_z).")
    print("\nRadial component S_r:")
    sympy.pprint(S_r_out, use_unicode=True)
    print("\nAxial component S_z:")
    sympy.pprint(S_z_out, use_unicode=True)
    print("-" * 40)

if __name__ == '__main__':
    compute_poynting_vector()