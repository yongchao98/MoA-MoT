import sympy

def solve_yukawa_ratio():
    """
    Calculates the ratio R for the given Yukawa theory.

    The calculation proceeds by determining the one-loop counter-terms
    delta_Z_x, delta_Z_m_x, and delta_Z_g in the MS-bar scheme.

    Let C be the common divergent factor C = g^2 / (32 * pi^2 * epsilon).
    All counter-terms will be expressed as multiples of C.
    """

    # 1. Fermion Self-Energy Calculation
    # The one-loop fermion self-energy Sigma(p) has a divergent part:
    # Sigma_div(p) = [g^2 / (16*pi^2*epsilon)] * (1/2 * p_slash + M_x)
    # Sigma_div(p) = (2*C) * (1/2 * p_slash + M_x) = C*p_slash + 2*C*M_x
    # The counter-term part of the Lagrangian is (dZ_x * p_slash - (dZ_x + dZ_m_x) * M_x).
    # To cancel the divergence, the counter-terms must equal Sigma_div.
    # dZ_x * p_slash = C * p_slash  => dZ_x = C
    dZ_x_coeff = 1
    # (dZ_x + dZ_m_x) * M_x = 2*C*M_x => dZ_x + dZ_m_x = 2*C
    # dZ_m_x = 2*C - dZ_x = 2*C - 1*C = 1*C
    dZ_m_x_coeff = 1

    # 2. Vertex Correction Calculation
    # The divergent part of the one-loop vertex correction, Lambda_div, is:
    # Lambda_div = -g^2 / (16*pi^2*epsilon) = -2*C
    # The vertex counter-term dZ_1 is given by dZ_1 = -Lambda_div.
    # dZ_1 = -(-2*C) = 2*C
    dZ_1_coeff = 2

    # 3. Yukawa Coupling Renormalization
    # The counter-term for the coupling, dZ_g, is related to the others by:
    # dZ_g = dZ_1 - dZ_x - (1/2)*dZ_phi
    # The problem provides the condition that dZ_phi = 0.
    dZ_phi_coeff = 0
    dZ_g_coeff = dZ_1_coeff - dZ_x_coeff - 0.5 * dZ_phi_coeff

    # 4. Calculate the final ratio R
    # The common factor C will cancel out in the ratio.
    R = dZ_x_coeff / (dZ_g_coeff + dZ_m_x_coeff)

    # Print the step-by-step derivation
    print("This script calculates the ratio R = delta_Z_x / (delta_Z_g + delta_Z_m_x)")
    print("Let C be the common factor g^2 / (32 * pi^2 * epsilon_uv).\n")

    print(f"1. From the fermion self-energy, we find the counter-term coefficients:")
    print(f"   delta_Z_x   = {dZ_x_coeff} * C")
    print(f"   delta_Z_m_x = {dZ_m_x_coeff} * C\n")

    print(f"2. From the vertex correction, we find the vertex counter-term coefficient:")
    print(f"   delta_Z_1 = {dZ_1_coeff} * C\n")

    print(f"3. Using the relation delta_Z_g = delta_Z_1 - delta_Z_x - 0.5 * delta_Z_phi")
    print(f"   and the given condition delta_Z_phi = 0:")
    print(f"   delta_Z_g = ({dZ_1_coeff}*C) - ({dZ_x_coeff}*C) - (0.5 * {dZ_phi_coeff}*C) = {dZ_g_coeff} * C\n")

    print("4. Now we compute the ratio R:")
    print(f"   R = (delta_Z_x) / (delta_Z_g + delta_Z_m_x)")
    print(f"   R = ({dZ_x_coeff} * C) / (({dZ_g_coeff} * C) + ({dZ_m_x_coeff} * C))")
    # Simplify by showing the coefficients
    print(f"   R = {dZ_x_coeff} / ({dZ_g_coeff} + {dZ_m_x_coeff})")
    print(f"   R = {dZ_x_coeff} / {dZ_g_coeff + dZ_m_x_coeff}")
    print(f"   R = {R}\n")

    print(f"The final calculated value for the ratio R is {R}.")

solve_yukawa_ratio()