def calculate_counter_term_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The ratio is defined as R = δZ_x / (δZ_g + δZ_{m_x}).
    The coefficients for the counter-terms at one-loop in the MS-bar scheme are:
    δZ_x   ∝ -1
    δZ_{m_x} ∝ +3
    δZ_g   ∝ -5
    """

    # Coefficients of the common factor g^2 / (32*pi^2*epsilon_bar)
    dzx_coeff = -1
    dzmx_coeff = 3
    dzg_coeff = -5

    # Calculate the numerator and denominator of the ratio R
    numerator = dzx_coeff
    denominator = dzg_coeff + dzmx_coeff

    # Calculate the final result
    if denominator == 0:
        result_str = "undefined (division by zero)"
    else:
        result = numerator / denominator
        result_str = f"{result}"

    # Print the calculation step-by-step
    print("The ratio R is calculated from the coefficients of the one-loop counter-terms.")
    print("Let the common factor be C = g^2 / (32*pi^2*epsilon_bar).")
    print(f"δZ_x = ({dzx_coeff}) * C")
    print(f"δZ_{m_x} = ({dzmx_coeff}) * C")
    print(f"δZ_g = ({dzg_coeff}) * C")
    print("\nNow, we compute the ratio R = δZ_x / (δZ_g + δZ_{m_x}):")
    print(f"R = ({dzx_coeff} * C) / (({dzg_coeff} * C) + ({dzmx_coeff} * C))")
    print("The common factor C cancels out, leaving the ratio of the coefficients:")
    print(f"R = ({numerator}) / (({dzg_coeff}) + ({dzmx_coeff}))")
    print(f"R = ({numerator}) / ({denominator})")
    print(f"R = {result_str}")

calculate_counter_term_ratio()