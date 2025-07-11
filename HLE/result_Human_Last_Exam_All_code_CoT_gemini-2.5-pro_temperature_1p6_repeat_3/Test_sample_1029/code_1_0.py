def compute_poynting_vector():
    """
    This function computes and displays the components of the Poynting vector S
    for a moving, charged cylindrical rod in an external electric field.

    The vector is expressed in cylindrical coordinates (r, phi, z).
    S = S_r * r_hat + S_z * z_hat

    Variables:
    E: Magnitude of the external uniform electric field (along the z-axis)
    rho: Uniform volume charge density of the rod
    v: Speed of the rod (along the z-axis)
    r: Radial distance from the center of the rod (r < R)
    epsilon_0: Permittivity of free space
    """

    # Construct the components of the Poynting vector as strings
    s_radial_component = "-(E * rho * v * r) / 2"
    s_axial_component = "(rho**2 * v * r**2) / (4 * epsilon_0)"

    print("The Poynting vector S is given by S = S_r * r_hat + S_z * z_hat, where:")
    print("-" * 50)

    # Print the radial component
    print(f"The radial component (S_r) is:")
    print(f"S_r = {s_radial_component}\n")
    # Outputting numbers in the equation
    print(f"Here, the number in the equation is: 2")


    # Print the axial component
    print(f"The axial component (S_z) is:")
    print(f"S_z = {s_axial_component}\n")
    # Outputting numbers in the equation
    print(f"Here, the number in the equation is: 4")

    # Print the full vector
    print("-" * 50)
    print("The full Poynting vector S is:")
    print(f"S = ({s_radial_component}) r_hat + ({s_axial_component}) z_hat")


if __name__ == '__main__':
    compute_poynting_vector()
