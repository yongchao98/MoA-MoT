def calculate_counter_term_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The ratio is defined as R = dZ_x / (dZ_g + dZ_m_x), where dZ_i are the
    one-loop counter-term coefficients in the MS-bar scheme.

    The calculation proceeds by determining the coefficients of the divergent
    term K = g^2 / (32*pi^2*epsilon) for each counter-term.
    """

    print("Step-by-step calculation of the ratio R = delta_Z_x / (delta_Z_g + delta_Z_m_x)")
    print("------------------------------------------------------------------------------------")
    print("Let K be the common factor for divergence: K = g^2 / (32 * pi^2 * epsilon)\n")

    # 1. Coefficient for the fermion field counter-term, delta_Z_x
    # From fermion self-energy, delta_Z_x = (1) * K
    dZ_x_coeff = 1.0
    print(f"1. The fermion field counter-term is delta_Z_x = {dZ_x_coeff} * K")

    # 2. Coefficient for the fermion mass counter-term, delta_Z_m_x
    # From fermion self-energy, delta_Z_m_x = (-3) * K
    dZ_m_coeff = -3.0
    print(f"2. The fermion mass counter-term is delta_Z_m_x = {dZ_m_coeff} * K")

    # 3. Coefficient for the vertex counter-term, delta_1
    # From the one-loop vertex correction diagram, delta_1 = 2 * K
    dZ1_coeff = 2.0
    print(f"3. The vertex correction counter-term is delta_1 = {dZ1_coeff} * K")

    # 4. Calculate the coefficient for the Yukawa coupling counter-term, delta_Z_g
    # Using the relation delta_1 = delta_Z_g + delta_Z_x (since delta_Z_phi = 0)
    dZ_g_coeff = dZ1_coeff - dZ_x_coeff
    print(f"4. The Yukawa coupling counter-term is delta_Z_g = delta_1 - delta_Z_x = ({dZ1_coeff} - {dZ_x_coeff}) * K = {dZ_g_coeff} * K\n")

    # 5. Compute the final ratio R
    print("5. Now, we compute the ratio R:")
    numerator = dZ_x_coeff
    denominator = dZ_g_coeff + dZ_m_coeff

    print(f"   R = (delta_Z_x) / (delta_Z_g + delta_Z_m_x)")
    print(f"   R = ({dZ_x_coeff} * K) / (({dZ_g_coeff} * K) + ({dZ_m_coeff} * K))")
    print("   The factor K cancels out.")
    print(f"   R = {numerator} / ({dZ_g_coeff} + {dZ_m_coeff})")
    print(f"   R = {numerator} / ({denominator})")

    # Final result
    R = numerator / denominator
    print(f"   R = {R}")

    return R

if __name__ == '__main__':
    final_ratio = calculate_counter_term_ratio()
    print("\n<<<{}>>>".format(final_ratio))