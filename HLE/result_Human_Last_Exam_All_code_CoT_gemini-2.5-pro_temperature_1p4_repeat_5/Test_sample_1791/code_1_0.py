import sympy

def derive_effective_interaction():
    """
    This function symbolically derives and prints the effective electron-electron
    interaction Hamiltonian mediated by phonons for a given wavevector q.
    """
    # Define the symbolic variables for the physical quantities.
    g = sympy.Symbol('g')          # Electron-phonon coupling constant
    q_j = sympy.Symbol('q_j')      # j-th component of the wavevector q
    m = sympy.Symbol('m')          # Mass of the ion
    w_q = sympy.Symbol('w_q')      # Phonon frequency for wavevector q
    n_q = sympy.Symbol('n_q')      # Electron density operator for wavevector q
    n_mq = sympy.Symbol('n_{-q}')  # Electron density operator for wavevector -q
    j = sympy.Symbol('j')          # Summation index for polarization/dimension
    d = sympy.Symbol('d')          # Number of polarizations/dimensions

    # --- Derivation Steps ---
    # The effective interaction H_eff(q) = V_eff(q) * n_q * n_{-q}, where V_eff(q) is the
    # effective potential derived by integrating out the phonons.
    # V_eff(q) is found from second-order perturbation theory or path integrals.
    # V_eff(q, omega) = -Sum_j |M_qj|^2 * D_ph(q, j, omega)
    # where M_qj is the interaction vertex and D_ph is the phonon propagator.

    # 1. From the interaction Hamiltonian, the vertex factor |M_qj|^2 is derived.
    # H_el-ph = g * Sum_j (i*q_j / sqrt(2*m*w_q)) * n_q * (a_{q,j} + a_{-q,j}^dagger)
    # The term in the parenthesis is the vertex M_qj.
    # Its absolute square is |M_qj|^2 = g^2 * q_j^2 / (2*m*w_q).
    abs_M_qj_sq = (g**2 * q_j**2) / (2 * m * w_q)

    # 2. The propagator for the phonon operator (a_{q,j} + a_{-q,j}^dagger) for zero frequency
    # transfer (static interaction) is D_ph(q, j, omega=0) = 2/w_q.
    D_ph_static = 2 / w_q

    # 3. The effective potential for a single polarization j is V_eff(q,j) = -|M_qj|^2 * D_ph.
    V_eff_qj = -abs_M_qj_sq * D_ph_static

    # --- Final Equation ---
    # The final Hamiltonian is constructed by summing the potential over all polarizations j
    # and multiplying by the electron density operators.
    
    print("The effective electron-electron interaction Hamiltonian for a given wavevector q is:")
    
    # We explicitly construct the string for the final equation to show each component clearly.
    # Each component of the formula is represented as a string.
    
    term_g_squared = "g^2"
    term_q_j_squared = "q_j^2"
    term_m = "m"
    term_w_q_squared = "w_q^2"
    term_n_q = "n_q"
    term_n_mq = "n_{-q}"

    # The potential term for a single polarization j.
    potential_term_string = f"- ( {term_g_squared} * {term_q_j_squared} ) / ( {term_m} * {term_w_q_squared} )"
    
    # The total potential is the sum over j.
    total_potential_string = f"Sum_j [ {potential_term_string} ]"
    
    # The final Hamiltonian expression.
    final_equation = f"H_eff(q) = ( {total_potential_string} ) * {term_n_q} * {term_n_mq}"

    print(final_equation)

if __name__ == '__main__':
    derive_effective_interaction()