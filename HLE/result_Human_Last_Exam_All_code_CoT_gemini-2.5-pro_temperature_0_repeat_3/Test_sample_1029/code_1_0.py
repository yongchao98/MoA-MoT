def solve_poynting_vector():
    """
    This function calculates and displays the formula for the Poynting vector
    for the given physics problem.

    The problem involves a long insulating cylindrical rod with:
    - Radius: R
    - Uniform volume charge density: rho
    - Axial speed: v
    - An external uniform electric field parallel to its axis: E

    The Poynting vector S is calculated in cylindrical coordinates (r, phi, z)
    with the z-axis along the rod's axis.
    """

    # The Poynting vector S is given by S = (1/mu_0) * (E_total x B)
    # We derive the expressions for E_total and B first.
    # E_total = E_external + E_rod = E*z_hat + E_r(r)*r_hat
    # B = B_phi(r)*phi_hat

    # After performing the cross product, S has two components: S_r and S_z.
    # S = S_r * r_hat + S_z * z_hat
    # The components are different inside (r < R) and outside (r > R) the rod.

    # Expressions for the Poynting vector components inside the rod (r < R)
    # The numbers 2 and 4 are the coefficients in the final equation.
    s_r_in = "-(E * rho * v * r) / 2"
    s_z_in = "(rho**2 * v * r**2) / (4 * epsilon_0)"

    # Expressions for the Poynting vector components outside the rod (r > R)
    # The numbers 2 and 4 are the coefficients in the final equation.
    s_r_out = "-(E * rho * v * R**2) / (2 * r)"
    s_z_out = "(rho**2 * v * R**4) / (4 * epsilon_0 * r**2)"

    # Print the final results in a structured format
    print("The Poynting vector S is computed for two regions:")
    print("The vector is expressed in cylindrical coordinates as S = S_r * r_hat + S_z * z_hat.")
    print("Where:")
    print("  E = magnitude of the external electric field")
    print("  rho = volume charge density of the rod")
    print("  v = speed of the rod")
    print("  R = radius of the rod")
    print("  r = radial distance from the rod's axis")
    print("  epsilon_0 = permittivity of free space")
    print("-" * 60)

    print("1. Inside the rod (for r < R):")
    print(f"   S_r = {s_r_in}")
    print(f"   S_z = {s_z_in}")
    print(f"   S_inside = ({s_r_in}) r_hat + ({s_z_in}) z_hat")

    print("-" * 60)

    print("2. Outside the rod (for r > R):")
    print(f"   S_r = {s_r_out}")
    print(f"   S_z = {s_z_out}")
    print(f"   S_outside = ({s_r_out}) r_hat + ({s_z_out}) z_hat")

# Execute the function to print the solution
solve_poynting_vector()