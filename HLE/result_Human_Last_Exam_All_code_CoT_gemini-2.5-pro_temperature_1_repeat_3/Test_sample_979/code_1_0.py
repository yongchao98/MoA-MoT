def print_magnetic_field_solution():
    """
    Prints the final derived expression for the magnetic field H(r, θ)
    in the two regions, inside and outside the spherical shell.
    """
    
    # Unicode characters for better readability
    mu_0 = "\u03BC\u2080"
    mu = "\u03BC"
    theta = "\u03B8"
    r_hat = "\u0072\u0302"
    theta_hat = "\u03B8\u0302"
    z_hat = "\u007A\u0302"

    # Expression for the inside field (0 < r < R)
    inside_field = f"H_in = ( (2 * {mu_0} / {mu}) * K\u2080 / (1 + (2 * {mu_0} / {mu})) ) * {z_hat}"

    # Expression for the outside field (r > R)
    outside_field = f"H_out = ( K\u2080 * R\u00b3 / ( (1 + (2 * {mu_0} / {mu})) * r\u00b3 ) ) * (2 * cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"

    print("The derived magnetic field H(r, \u03B8) is given by a piecewise function:")
    print("-" * 60)

    print("\nFor the region inside the sphere (0 < r < R):")
    print(inside_field)

    print("\nFor the region outside the sphere (r > R):")
    print(outside_field)
    print("-" * 60)
    
    print("\nThis corresponds to Answer Choice E.")

if __name__ == "__main__":
    print_magnetic_field_solution()