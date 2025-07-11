def poynting_vector_calculation():
    """
    This function prints the symbolic expression for the Poynting vector
    for a moving, charged cylindrical rod in an external electric field.

    The derivation involves finding the electric and magnetic fields first,
    then computing their cross product.

    Variables:
    E: Magnitude of the uniform external electric field (along z-axis)
    rho: Uniform volume charge density of the rod
    v: Speed of the rod (along z-axis)
    R: Radius of the rod
    r: Radial distance from the axis of the rod
    epsilon_0: Permittivity of free space
    mu_0: Permeability of free space

    Unit vectors:
    r_hat: Radial unit vector
    z_hat: Axial unit vector
    """

    # Components of the Poynting vector S = S_r * r_hat + S_z * z_hat

    # For r < R (inside the rod)
    # E_total = (rho*r / (2*epsilon_0)) * r_hat + E * z_hat
    # B = (mu_0*rho*v*r / 2) * phi_hat
    # S = (1/mu_0) * (E_total x B)
    s_r_inside = "-(E * rho * v * r / 2)"
    s_z_inside = "(rho**2 * v * r**2 / (4 * epsilon_0))"

    # For r > R (outside the rod)
    # E_total = (rho*R**2 / (2*epsilon_0*r)) * r_hat + E * z_hat
    # B = (mu_0*rho*v*R**2 / (2*r)) * phi_hat
    # S = (1/mu_0) * (E_total x B)
    s_r_outside = "-(E * rho * v * R**2 / (2 * r))"
    s_z_outside = "(rho**2 * v * R**4 / (4 * epsilon_0 * r**2))"

    print("The Poynting vector S describes the energy flux density and is given by S = S_r * r_hat + S_z * z_hat.")
    print("The vector components depend on the radial distance 'r' from the axis.\n")

    print("For r < R (inside the rod):")
    print(f"S_r = {s_r_inside}")
    print(f"S_z = {s_z_inside}")
    print(f"So, the Poynting vector is: S = {s_r_inside} r_hat + {s_z_inside} z_hat\n")

    print("For r > R (outside the rod):")
    print(f"S_r = {s_r_outside}")
    print(f"S_z = {s_z_outside}")
    print(f"So, the Poynting vector is: S = {s_r_outside} r_hat + {s_z_outside} z_hat")

if __name__ == '__main__':
    poynting_vector_calculation()