def calculate_counterterm_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The ratio is defined as R = dZx / (dZg + dZmx), where:
    - dZx is the fermion field renormalization counter-term.
    - dZg is the Yukawa coupling renormalization counter-term.
    - dZmx is the fermion mass renormalization counter-term.

    In the MS-bar scheme, with the condition dZ_phi = 0, we find:
    - dZx = C
    - dZg = -C
    - dZmx = -3C
    where C is the common factor g^2 / (32 * pi^2 * epsilon).
    """

    # Define a symbolic representation for the common factor
    C_str = "g^2 / (32 * pi^2 * epsilon)"

    # Define the coefficients of the counter-terms in units of C
    dZx_coeff = 1.0
    dZg_coeff = -1.0
    dZmx_coeff = -3.0

    print("Step 1: Expressing the one-loop counter-terms.")
    print(f"Let C = {C_str}")
    print(f"From one-loop calculations in the MS-bar scheme:")
    print(f"delta Z_x   = ({dZx_coeff}) * C")
    print(f"delta Z_g   = ({dZg_coeff}) * C")
    print(f"delta Z_m_x = ({dZmx_coeff}) * C")
    print("-" * 30)

    # Calculate the denominator
    denominator_coeff = dZg_coeff + dZmx_coeff
    
    # Calculate the ratio R
    R = dZx_coeff / denominator_coeff

    print("Step 2: Calculating the ratio R = delta Z_x / (delta Z_g + delta Z_m_x).")
    print("The final equation for the ratio is:")
    print(f"R = ({dZx_coeff} * C) / (({dZg_coeff} * C) + ({dZmx_coeff} * C))")
    print(f"R = ({dZx_coeff} * C) / ({denominator_coeff} * C)")
    print(f"R = {dZx_coeff} / {denominator_coeff}")
    print("-" * 30)
    
    print("Step 3: Final Result.")
    print(f"The calculated value of the ratio R is: {R}")

if __name__ == "__main__":
    calculate_counterterm_ratio()