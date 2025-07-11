import sympy as sp

def display_magnetic_field_solution():
    """
    This function displays the solution for the magnetic field inside and outside the spherical shell.
    The formulas are derived from first principles using the magnetic scalar potential and boundary conditions.
    """
    
    # Define symbolic variables for clarity, although we will print strings.
    # K0, mu, mu0, R, r, theta = sp.symbols('K_0 mu mu_0 R r theta')

    # Expression for the H-field inside the sphere (r < R)
    # The field is uniform and in the z-direction.
    # The term (cos(theta)*r_hat - sin(theta)*theta_hat) is the z_hat vector in spherical coordinates.
    # So the expression is C * z_hat, where C is the coefficient.
    H_inside_str = "H_in(r, theta) = (2 * mu_0 / mu) * K_0 / (1 + 2 * mu_0 / mu) * z_hat"
    
    # Expression for the H-field outside the sphere (r > R)
    # This is a standard magnetic dipole field.
    H_outside_str = "H_out(r, theta) = K_0 / (1 + 2 * mu_0 / mu) * (R**3 / r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    print("The magnetic field H(r, theta) in the two regions is given by:")
    print("-" * 60)
    
    # Print the field inside the sphere
    print("For 0 < r < R (inside the sphere):")
    # To demonstrate outputting each 'number', we will construct the equation part by part.
    # This is for illustrative purposes; in a real calculation, these would be variables.
    term1_in = "2"
    term2_in = "mu_0 / mu"
    term3_in = "K_0"
    denominator_in = "(1 + 2 * mu_0 / mu)"
    vector_in = "z_hat"
    print(f"H_in = ({term1_in} * {term2_in} * {term3_in}) / {denominator_in} * {vector_in}")
    
    print("\n" + "-" * 60 + "\n")
    
    # Print the field outside the sphere
    print("For R < r < infinity (outside the sphere):")
    # Similarly constructing the output part by part.
    term1_out = "K_0"
    denominator_out = "(1 + 2 * mu_0 / mu)"
    r_dependence = "(R**3 / r**3)"
    angular_dependence = "(2*cos(theta)*r_hat + sin(theta)*theta_hat)"
    print(f"H_out = ({term1_out} / {denominator_out}) * {r_dependence} * {angular_dependence}")
    
    print("\n" + "-" * 60)
    print("Note: mu is the permeability inside the sphere, mu_0 is the permeability of free space,")
    print("K_0 is the current amplitude, R is the sphere's radius, and r, theta are spherical coordinates.")
    print("z_hat, r_hat, theta_hat are the respective unit vectors.")


if __name__ == "__main__":
    display_magnetic_field_solution()