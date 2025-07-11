def print_magnetic_field_solution():
    """
    This function prints the derived equations for the magnetic field H
    in the two specified regions.
    """

    # Print the field equation for the region inside the shield (0 < r < Rp)
    print("In the region 0 < r < R_p:")
    # Using 'i_r' and 'i_theta' for the unit vectors
    h_field_inner = (
        "  H = M_0 * ((2 * R_p**3 + R**3) / (3 * R**3)) * "
        "(-cos(theta) * i_r + sin(theta) * i_theta)"
    )
    print(h_field_inner)
    print("\n")  # Add a newline for better separation

    # Print the field equation for the region between the shield and the conductor (Rp < r < R)
    print("In the region R_p < r < R:")
    h_field_outer_r_component = (
        "  H_r = (-2 * M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r"
    )
    h_field_outer_theta_component = (
        "  H_theta = (M_0 / 3) * [ 2 * (R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"
    )

    # Combining both components for the total field H
    print("  The components of the magnetic field H are:")
    print(h_field_outer_r_component)
    print(h_field_outer_theta_component)
    print("  So, H = H_r + H_theta")

if __name__ == "__main__":
    print_magnetic_field_solution()
