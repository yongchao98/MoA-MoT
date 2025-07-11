def compute_poynting_vector():
    """
    This function calculates and prints the symbolic expression for the Poynting vector
    both inside and outside the described cylindrical rod.
    """
    # --- Variable Definitions ---
    # E: Magnitude of the external electric field (along z-axis)
    # ρ: Uniform volume charge density of the rod
    # v: Speed of the rod (along z-axis)
    # R: Radius of the rod
    # s: Radial distance from the axis of the rod
    # ε₀: Permittivity of free space
    # The unit vectors are s_hat (radial) and z_hat (axial).

    print("The Poynting vector S is calculated for two regions: inside the rod (s < R) and outside the rod (s > R).")
    print("The vector has components in the radial (s_hat) and axial (z_hat) directions.\n")

    # --- Inside the rod (s < R) ---
    print("--- Inside the rod (s < R) ---")
    
    # Derivations give:
    # E_total = (ρ*s / (2*ε₀)) s_hat + E z_hat
    # B = (μ₀*ρ*v*s / 2) φ_hat
    # S = (1/μ₀) * (E_total x B)
    # S = - (E*ρ*v*s / 2) s_hat + (ρ²*v*s² / (4*ε₀)) z_hat
    
    s_radial_in_expr = "-(E * ρ * v * s / 2)"
    s_axial_in_expr = "(ρ**2 * v * s**2) / (4 * ε₀)"
    
    print("The Poynting vector S_inside is given by:")
    print(f"S_inside = {s_radial_in_expr} s_hat + {s_axial_in_expr} z_hat\n")
    print("Each term in the final equation is:")
    print(f"  Radial Component: {s_radial_in_expr}")
    print(f"  Axial Component:  {s_axial_in_expr}")
    print("-" * 50)

    # --- Outside the rod (s > R) ---
    print("\n--- Outside the rod (s > R) ---")
    
    # Derivations give:
    # E_total = (ρ*R² / (2*ε₀*s)) s_hat + E z_hat
    # B = (μ₀*ρ*v*R² / (2*s)) φ_hat
    # S = (1/μ₀) * (E_total x B)
    # S = - (E*ρ*v*R² / (2*s)) s_hat + (ρ²*v*R⁴ / (4*ε₀*s²)) z_hat
    
    s_radial_out_expr = "-(E * ρ * v * R**2) / (2 * s)"
    s_axial_out_expr = "(ρ**2 * v * R**4) / (4 * ε₀ * s**2)"
    
    print("The Poynting vector S_outside is given by:")
    print(f"S_outside = {s_radial_out_expr} s_hat + {s_axial_out_expr} z_hat\n")
    print("Each term in the final equation is:")
    print(f"  Radial Component: {s_radial_out_expr}")
    print(f"  Axial Component:  {s_axial_out_expr}")

# Execute the function to print the solution
compute_poynting_vector()