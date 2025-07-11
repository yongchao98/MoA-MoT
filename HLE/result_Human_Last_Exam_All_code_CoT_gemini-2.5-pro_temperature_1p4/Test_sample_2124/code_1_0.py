def calculate_yukawa_ratio():
    """
    Calculates the ratio R = δZ_x / (δZ_g + δZ_m_x) for a Yukawa theory.

    The calculation is based on one-loop counter-terms in the MS-bar scheme,
    under the condition that the scalar field renormalization δZ_φ is zero.

    We express the counter-terms in units of a common factor C = g^2 / (16*pi^2*epsilon).
    - δZ_x = -1/2 * C
    - δZ_m_x = -1 * C
    - δZ_g = 0 (due to δZ_1 = δZ_x and the given condition δZ_φ = 0)
    """

    # Coefficients of the common factor C = g^2/(16*pi^2*epsilon)
    dz_x_coeff = -0.5
    dz_mx_coeff = -1.0
    dz_g_coeff = 0.0

    # The ratio R is independent of the common factor C, so we can use the coefficients directly.
    numerator = dz_x_coeff
    denominator = dz_g_coeff + dz_mx_coeff
    
    # Calculate the ratio
    if denominator == 0:
        R = float('inf') if numerator != 0 else float('nan')
    else:
        R = numerator / denominator

    # Print the details of the calculation
    print("This script calculates the ratio R = δZ_x / (δZ_g + δZ_{m_x}).")
    print("The counter-terms can be expressed in terms of a common factor C = g^2 / (16*pi^2*epsilon).")
    print(f"The coefficient for δZ_x is: {dz_x_coeff}")
    print(f"The coefficient for δZ_g is: {dz_g_coeff}")
    print(f"The coefficient for δZ_{m_x} is: {dz_mx_coeff}")
    print("\nThe final equation for R using these coefficients is:")
    print(f"R = ({dz_x_coeff}) / ({dz_g_coeff} + ({dz_mx_coeff}))")
    print(f"R = {dz_x_coeff} / {denominator}")
    print(f"R = {R}")
    print("\n<<<{}>>>".format(R))


if __name__ == "__main__":
    calculate_yukawa_ratio()