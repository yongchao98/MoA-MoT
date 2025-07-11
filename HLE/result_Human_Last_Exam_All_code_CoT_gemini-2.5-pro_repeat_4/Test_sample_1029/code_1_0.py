def compute_poynting_vector():
    """
    This function derives and prints the symbolic expression for the Poynting vector
    for the given physics problem.
    
    The final expression is split into two cases: inside the rod (r < R) and
    outside the rod (r > R).
    
    The equation is S = (1/mu_0) * (E_total x B)
    E_total = E_external + E_rod
    
    The code constructs the final string expressions and prints them.
    """
    
    # Symbolic representation of the variables
    # E: External electric field magnitude
    # rho: Volume charge density
    # v: Speed of the rod
    # R: Radius of the rod
    # r: Radial distance from the axis
    # epsilon_0: Permittivity of free space
    # r_hat, z_hat: Unit vectors
    
    # Poynting vector inside the rod (r < R)
    # The radial component has a coefficient of -1/2
    # The axial component has a coefficient of 1/4
    s_inside_radial_coeff_num = -1
    s_inside_radial_coeff_den = 2
    s_inside_axial_coeff_num = 1
    s_inside_axial_coeff_den = 4
    
    s_inside = (f"S = ({s_inside_radial_coeff_num}/{s_inside_radial_coeff_den}) * E * rho * v * r * r_hat  +  "
                f"({s_inside_axial_coeff_num}/{s_inside_axial_coeff_den}) * (rho^2 * v * r^2 / epsilon_0) * z_hat")
    
    # Poynting vector outside the rod (r > R)
    # The radial component has a coefficient of -1/2
    # The axial component has a coefficient of 1/4
    s_outside_radial_coeff_num = -1
    s_outside_radial_coeff_den = 2
    s_outside_axial_coeff_num = 1
    s_outside_axial_coeff_den = 4
    
    s_outside = (f"S = ({s_outside_radial_coeff_num}/{s_outside_radial_coeff_den}) * (E * rho * v * R^2 / r) * r_hat  +  "
                 f"({s_outside_axial_coeff_num}/{s_outside_axial_coeff_den}) * (rho^2 * v * R^4 / (epsilon_0 * r^2)) * z_hat")
                 
    print("The Poynting vector S has two components: a radial one (r_hat) and an axial one (z_hat).")
    print("It depends on whether the point is inside or outside the rod.")
    print("-" * 70)
    
    print("Inside the rod (for r < R):\n")
    print(s_inside)
    
    print("\n" + "-" * 70)
    
    print("Outside the rod (for r > R):\n")
    print(s_outside)
    print("-" * 70)

if __name__ == "__main__":
    compute_poynting_vector()