import sympy

def print_field_equation():
    """
    This function prints the derived field equation for the given theory of gravity.
    Each numerical coefficient in the equation is explicitly shown.
    """
    # Define symbols for the equation using LaTeX syntax
    sqrt_g = r"\sqrt{-g}"
    partial_alpha = r"\partial_{\alpha}"
    P_alpha_munu = r"P^\alpha{}_{\mu\nu}"
    P_mu_alpha_beta = r"P_{\mu\alpha\beta}"
    Q_nu_alpha_beta = r"Q_\nu{}^{\alpha\beta}"
    Q_alpha_beta_mu = r"Q^{\alpha\beta}{}_\mu"
    P_alpha_beta_nu = r"P_{\alpha\beta\nu}"
    Q_scalar = "Q"
    g_munu = r"g_{\mu\nu}"
    T_munu = r"T_{\mu\nu}"
    constants = r"\frac{8\pi G}{c^4}"

    # Build the equation string part by part
    
    # Term 1: Derivative of the superpotential
    term1_coeff = -2
    term1_str = f"\\frac{{{term1_coeff}}}{{{sqrt_g}}}{partial_alpha}({sqrt_g}{P_alpha_munu})"
    
    # Term 2: P*Q term
    term2_coeff = -2
    term2_str = f"{term2_coeff}{P_mu_alpha_beta} {Q_nu_alpha_beta}"
    
    # Term 3: Q*P term
    term3_coeff = 1
    term3_str = f"+{term3_coeff}{Q_alpha_beta_mu} {P_alpha_beta_nu}"
    
    # Term 4: Scalar term
    term4_coeff_num = -1
    term4_coeff_den = 2
    term4_str = f"\\frac{{{term4_coeff_num}}}{{{term4_coeff_den}}}{Q_scalar}{g_munu}"
    
    # Right-hand side (RHS)
    rhs_coeff = 8
    rhs_str = f"\\frac{{{rhs_coeff}\\pi G}}{{c^4}} {T_munu}"
    
    # Combine all parts into the final equation string
    final_equation = f"{term1_str} {term2_str} {term3_str} {term4_str} = {rhs_str}"
    
    print("The derived field equation is:")
    # For better display in environments that support LaTeX rendering like Jupyter notebooks
    try:
        from IPython.display import display, Math
        display(Math(final_equation))
    except ImportError:
        print(final_equation)

    print("\nNumerical coefficients in the final equation:")
    print(f"Coefficient of the derivative term: {term1_coeff}")
    print(f"Coefficient of the P*Q term: {term2_coeff}")
    print(f"Coefficient of the Q*P term: {term3_coeff}")
    print(f"Coefficient of the Qg term: {term4_coeff_num}/{term4_coeff_den}")
    print(f"Numerical constant on the RHS: {rhs_coeff}")

if __name__ == "__main__":
    print_field_equation()