def display_poynting_vector_solution():
    """
    This script prints the derived symbolic expression for the Poynting vector S.
    The final equation involves the following physical quantities and constants:

    Variables:
    E:         Magnitude of the external uniform electric field along the axis.
    rho:       Uniform volume charge density of the rod.
    v:         Speed of the rod along its axis.
    R:         Radius of the cylindrical rod.
    r:         Radial distance from the center of the rod.

    Constants:
    epsilon_0: Permittivity of free space.

    Vectors:
    r_hat:     Unit vector in the radial direction, away from the axis.
    z_hat:     Unit vector along the axis of the rod.
    """
    print("The Poynting vector S is calculated in two regions:")
    print("="*60)

    # Solution for inside the rod (r < R)
    print("1. Inside the rod (where r < R):")
    # Equation
    equation_inside = "S_vector = - (E * rho * v * r / 2) * r_hat   +   (rho^2 * v * r^2 / (4 * epsilon_0)) * z_hat"
    print(equation_inside)
    
    # Components breakdown
    print("\n   The radial component (S_r) is: -(E * rho * v * r / 2)")
    print("   The axial component (S_z) is:   (rho^2 * v * r^2 / (4 * epsilon_0))")

    print("\n" + "="*60 + "\n")

    # Solution for outside the rod (r > R)
    print("2. Outside the rod (where r > R):")
    # Equation
    equation_outside = "S_vector = - (E * rho * v * R^2 / (2 * r)) * r_hat   +   (rho^2 * v * R^4 / (4 * epsilon_0 * r^2)) * z_hat"
    print(equation_outside)

    # Components breakdown
    print("\n   The radial component (S_r) is: -(E * rho * v * R^2 / (2 * r))")
    print("   The axial component (S_z) is:   (rho^2 * v * R^4 / (4 * epsilon_0 * r^2))")
    print("="*60)

if __name__ == '__main__':
    display_poynting_vector_solution()