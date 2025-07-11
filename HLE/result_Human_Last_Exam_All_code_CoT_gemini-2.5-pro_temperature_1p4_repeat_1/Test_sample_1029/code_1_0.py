def display_poynting_vector():
    """
    This function prints the derived symbolic expression for the Poynting vector S
    both inside (r < R) and outside (r > R) the moving cylindrical rod.

    The Poynting vector S has two components:
    - S_r: The radial component, representing energy flow into/out of the rod's axis.
    - S_z: The axial component, representing energy flow along the rod's axis.
    There is no azimuthal component (S_phi = 0).

    Variables used in the equations:
    E: Magnitude of the external uniform electric field
    ρ: Uniform volume charge density of the rod
    v: Speed of the rod along its axis
    R: Radius of the rod
    r: Radial distance from the center of the rod
    ε₀: Permittivity of free space
    μ₀: Permeability of free space
    """

    # --- Inside the rod (r < R) ---
    print("--- Poynting Vector Inside the Rod (r < R) ---")
    
    # S = S_r * r_hat + S_z * z_hat
    Sr_inside = "- (E * ρ * v * r) / 2"
    Sz_inside = "+ (ρ**2 * v * r**2) / (4 * ε₀)"
    
    print(f"The Poynting vector S is composed of a radial component (S_r) and an axial component (S_z).")
    print(f"S_r (radial component) = {Sr_inside}")
    print(f"S_z (axial component) = {Sz_inside}\n")
    print(f"Full Vector: S = {Sr_inside} r_hat {Sz_inside} z_hat")
    
    print("\n" + "="*50 + "\n")

    # --- Outside the rod (r > R) ---
    print("--- Poynting Vector Outside the Rod (r > R) ---")
    
    # S = S_r * r_hat + S_z * z_hat
    Sr_outside = "- (E * ρ * v * R**2) / (2 * r)"
    Sz_outside = "+ (ρ**2 * v * R**4) / (4 * ε₀ * r**2)"
    
    print(f"The Poynting vector S is composed of a radial component (S_r) and an axial component (S_z).")
    print(f"S_r (radial component) = {Sr_outside}")
    print(f"S_z (axial component) = {Sz_outside}\n")
    print(f"Full Vector: S = {Sr_outside} r_hat {Sz_outside} z_hat")


# Execute the function to display the results
if __name__ == "__main__":
    display_poynting_vector()
