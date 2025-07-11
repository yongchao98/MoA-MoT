def solve_magnetic_field_sphere():
    """
    This function provides the final analytical solution for the magnetic field H(r, θ)
    inside and outside a spherical shell with a given surface current.
    The derivation involves solving Laplace's equation for the magnetic scalar potential
    and applying the appropriate boundary conditions at the surface of the sphere.
    """

    # --- Final Expression Components ---
    # The term '1' in the denominators
    one = 1

    # The coefficient '2' appearing in several places
    two = 2
    
    # Inside field H_in
    H_in_coeff = f"({two} * μ₀ / μ) * K₀ / ({one} + {two} * μ₀ / μ)"
    H_in_vector = "ẑ"
    H_in_expression = f"{H_in_coeff} {H_in_vector}"

    # Outside field H_out
    H_out_coeff = f"K₀ / ({one} + {two} * μ₀ / μ)"
    H_out_radial_dep = "R³ / r³"
    H_out_angular_dep = f"({two} * cos(θ) * r̂ + sin(θ) * θ̂)"
    H_out_expression = f"{H_out_coeff} * {H_out_radial_dep} * {H_out_angular_dep}"

    # --- Print the solution ---
    print("The final derived magnetic field H(r, θ) is:")
    print("H(r, θ) =")
    print("  {")
    print(f"  |  {H_in_expression}      for 0 < r < R")
    print(f"  |")
    print(f"  |  {H_out_expression}      for r > R")
    print("  }")
    print("\nThis result corresponds to answer choice E.")

if __name__ == '__main__':
    solve_magnetic_field_sphere()
