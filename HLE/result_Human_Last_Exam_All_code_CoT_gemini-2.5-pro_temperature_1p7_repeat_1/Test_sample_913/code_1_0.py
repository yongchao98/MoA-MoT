def display_final_equations():
    """
    This function prints the derived expressions for the electric field
    in the two specified regions. The expressions correspond to option B.
    """
    
    # Using unicode for better readability of greek letters and subscripts
    P_0 = "P\u2080"
    epsilon_0 = "\u03B5\u2080"
    R_p = "R\u209a"
    r_hat = "\u0155"
    theta_hat = "\u03B8\u0302"
    
    # Expression for the z-unit vector in spherical components
    z_hat_expr_str = f"(cos(\u03B8) {r_hat} - sin(\u03B8) {theta_hat})"
    
    # Expression for the dipole field vector part
    dipole_vector_str = f"(2*cos(\u03B8) {r_hat} + sin(\u03B8) {theta_hat})"

    print("The derived electric field expressions are as follows:")
    print("=" * 80)
    
    # Print the field for r < R_p
    print(f"For r < {R_p} (inside the sensor):")
    print(f"  \u20d7 = - ({P_0} / (3 * {epsilon_0})) * (1 - ({R_p}/R)\u00b3) * {z_hat_expr_str}")
    
    print("\n" + "=" * 80)
    
    # Print the field for R_p < r < R
    print(f"For {R_p} < r < R (in the free space):")
    print(f"  \u20d7 = ({P_0} / (3 * {epsilon_0})) * ({R_p}/R)\u00b3 * {z_hat_expr_str} + ({P_0} * {R_p}\u00b3 / (3 * {epsilon_0} * r\u00b3)) * {dipole_vector_str}")
    
    print("\n" + "=" * 80)
    print("These results match answer choice B.")

if __name__ == '__main__':
    display_final_equations()