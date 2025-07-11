def calculate_counter_term_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in a Yukawa theory.

    The calculation is based on the following one-loop results in the MS-bar scheme:
    - delta_Z_x = g^2 / (32 * pi^2 * epsilon)
    - delta_Z_m_x = -g^2 / (16 * pi^2 * epsilon)
    - delta_Z_g = g^2 / (32 * pi^2 * epsilon)

    The ratio R is defined as R = delta_Z_x / (delta_Z_g + delta_Z_m_x).
    The common factor C = g^2 / (32 * pi^2 * epsilon) will cancel out.
    We can perform the calculation using the numerical coefficients of C.
    """

    # Let C be the common symbolic factor g^2 / (32 * pi^2 * epsilon)
    # We will use its numerical coefficients for the calculation.
    
    # Coefficient for delta_Z_x
    dZx_coeff = 1.0
    
    # Coefficient for delta_Z_m_x. Note that 1/(16*pi^2) = 2 * 1/(32*pi^2)
    dZmx_coeff = -2.0
    
    # Coefficient for delta_Z_g
    dZg_coeff = 1.0

    print("This script calculates the ratio R = delta_Z_x / (delta_Z_g + delta_Z_m_x)")
    print("The counter-terms are expressed in units of C = g^2 / (32 * pi^2 * epsilon).")
    print("-" * 30)
    
    # Output each number in the final equation
    print(f"Value of delta_Z_x = {dZx_coeff} * C")
    print(f"Value of delta_Z_g = {dZg_coeff} * C")
    print(f"Value of delta_Z_m_x = {dZmx_coeff} * C")
    print("-" * 30)

    # Calculate the numerator and denominator of R
    numerator = dZx_coeff
    denominator = dZg_coeff + dZmx_coeff
    
    # Calculate the final ratio R
    R = numerator / denominator

    print(f"The equation for R is: R = {numerator} / ({dZg_coeff} + ({dZmx_coeff}))")
    print(f"R = {numerator} / {denominator}")
    print(f"The final result for the ratio R is: {R}")

calculate_counter_term_ratio()
<<< -1.0 >>>