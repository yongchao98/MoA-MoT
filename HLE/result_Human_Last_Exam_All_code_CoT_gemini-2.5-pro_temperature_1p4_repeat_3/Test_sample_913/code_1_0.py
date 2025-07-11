import sympy as sp

def display_electric_field_solution():
    """
    This function prints the derived electric field expressions for the given
    electrostatics problem.
    """

    # Define symbolic representations for readability
    P0 = sp.Symbol('P_0')
    eps0 = sp.Symbol('varepsilon_0')
    Rp = sp.Symbol('R_p')
    R = sp.Symbol('R')
    r = sp.Symbol('r')
    theta = sp.Symbol('theta')
    r_hat = sp.Symbol('r_hat')
    theta_hat = sp.Symbol('theta_hat')

    # Expression for the z-hat unit vector in spherical components
    z_hat_expr = f"(cos({theta}) * {r_hat} - sin({theta}) * {theta_hat})"

    # Expression for the dipole field term in spherical components
    dipole_field_expr = f"(2*cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"

    # Electric field inside the sensor (r < R_p)
    # E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * z_hat
    E_in_str = (
        f"E_vec = - ( {P0} / (3 * {eps0}) ) * (1 - ({Rp}/{R})**3) * "
        f"{z_hat_expr}"
    )

    # Electric field in the free space (R_p < r < R)
    # E_out = (P_0 / (3*eps0)) * (R_p/R)**3 * z_hat + (P_0*R_p**3 / (3*eps0*r**3)) * (dipole_term)
    E_out_str = (
        f"E_vec = ( {P0} / (3 * {eps0}) ) * ({Rp}/{R})**3 * {z_hat_expr}"
        f" + ( {P0} * {Rp}**3 / (3 * {eps0} * {r}**3) ) * {dipole_field_expr}"
    )

    print("The final expressions for the electric field are:")
    print("="*60)
    print("For r < R_p (inside the sensor):")
    print(E_in_str)
    print("\n" + "="*60)
    print("For R_p < r < R (in the free space):")
    print(E_out_str)
    print("="*60)
    print("\nThese expressions correspond to Answer Choice B.")

if __name__ == '__main__':
    display_electric_field_solution()
