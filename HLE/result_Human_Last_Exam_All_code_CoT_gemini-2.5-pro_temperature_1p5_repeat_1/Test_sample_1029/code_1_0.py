def display_poynting_vector_formulas():
    """
    This function prints the derived formulas for the Poynting vector (S)
    both inside and outside the described cylindrical rod.
    
    The formulas are expressed in terms of the given physical quantities:
    R: radius of the rod
    rho: volume charge density
    E: magnitude of the external electric field along the z-axis
    v: speed of the rod along the z-axis
    r: radial distance from the center of the rod
    epsilon_0: permittivity of free space
    
    The vector components are given in a cylindrical coordinate system (r_hat, z_hat).
    """

    # Formula for inside the rod (r < R)
    s_radial_in = "- (E * rho * v * r) / 2"
    s_axial_in = "+ (rho**2 * v * r**2) / (4 * epsilon_0)"
    
    print("Poynting Vector inside the rod (r < R):")
    print("S_in = ({}) r_hat {} z_hat".format(s_radial_in, s_axial_in))
    print("-" * 50)

    # Formula for outside the rod (r > R)
    s_radial_out = "- (E * rho * v * R**2) / (2 * r)"
    s_axial_out = "+ (rho**2 * v * R**4) / (4 * epsilon_0 * r**2)"

    print("Poynting Vector outside the rod (r > R):")
    print("S_out = ({}) r_hat {} z_hat".format(s_radial_out, s_axial_out))
    print("\n")
    print("Final answer in requested format:")
    # The final answer format is quite long, so we will print it on two lines.
    final_answer_inside = "For r < R, S = ({}) r_hat + ({}) z_hat".format(s_radial_in, s_axial_in.replace("+ ", ""))
    final_answer_outside = "For r > R, S = ({}) r_hat + ({}) z_hat".format(s_radial_out, s_axial_out.replace("+ ", ""))
    print(f"<<<{final_answer_inside}; {final_answer_outside}>>>")

if __name__ == "__main__":
    display_poynting_vector_formulas()
