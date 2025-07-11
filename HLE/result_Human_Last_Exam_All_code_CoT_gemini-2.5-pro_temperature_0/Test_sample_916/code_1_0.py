def print_force_equation():
    """
    Prints the components of the derived instantaneous force equation.
    """
    
    # Define the symbols used in the equation as strings
    term_force = "f_x(t)"
    term_equals = "="
    term_sign = "-"
    term_geometry = "2\pi R N"
    term_numerator_start = "("
    term_mu_temp = "\mu_0 (1 - \alpha_T (T - T_0))"
    term_mmf = "N_0 I_0"
    term_ac_current = "i_0 \sin(\omega t)"
    term_numerator_end = ")"
    term_denominator_start = "/"
    term_denominator = "g^2 (1 + \frac{\mu_0 N_0 I_0}{g B_s})"

    # Print the equation piece by piece
    print("The instantaneous force f_x(t) is given by the equation:")
    print(f"{term_force} {term_equals} {term_sign} {term_geometry} * [{term_numerator_start}{term_mu_temp} * {term_mmf} * {term_ac_current}{term_numerator_end} {term_denominator_start} {term_denominator}]")

print_force_equation()