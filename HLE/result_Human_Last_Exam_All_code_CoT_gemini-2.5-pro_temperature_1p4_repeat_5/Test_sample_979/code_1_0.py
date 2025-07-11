def solve_magnetic_field():
    """
    This function presents the derived expressions for the magnetic field H
    inside and outside the spherical shell.
    """
    
    # Define common terms as strings for clarity
    K0 = "K_0"
    mu = "mu"
    mu0 = "mu_0"
    R_cubed = "R^3"
    r_cubed = "r^3"

    # Common denominator factor
    denominator = f"(1 + 2*{mu0}/{mu})"

    # Expression for the magnetic field inside the sphere (0 < r < R)
    H_in = f"H_in(r, theta) = (2 * {mu0} / {mu}) * {K0} / {denominator} * z_hat"

    # Expression for the magnetic field outside the sphere (r > R)
    H_out = f"H_out(r, theta) = ({K0} / {denominator}) * ({R_cubed}/{r_cubed}) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    # Print the final combined solution
    print("The magnetic field H(r, theta) is given by:")
    print("{")
    print(f"  {H_in}  \t for 0 < r < R")
    print(f"  {H_out}  \t for r > R")
    print("}")
    print("\nThis corresponds to Answer Choice E.")

solve_magnetic_field()