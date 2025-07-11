import sympy

def derive_effective_interaction():
    """
    This function symbolically derives the effective electron-electron interaction
    mediated by phonons.
    """

    # --- 1. Define Symbolic Variables ---
    # We use sympy for symbolic representation, which makes the equations clear.
    g = sympy.Symbol('g')              # electron-phonon coupling constant
    m = sympy.Symbol('m')              # mass (from the coupling term)
    q_vec = sympy.Symbol('q')          # phonon wavevector
    q_j = sympy.Symbol('q_j')          # j-th component of the wavevector q
    w_q = sympy.Symbol('w_q')          # phonon frequency, omega_q
    omega_n = sympy.Symbol('omega_n')  # Matsubara frequency for bosons
    rho_q = sympy.Symbol('rho_q')      # electron density operator at wavevector q
    
    print("### Derivation of the Effective Electron-Electron Interaction ###\n")

    # --- 2. Starting Hamiltonians ---
    print("1. The starting Hamiltonians are:")
    print(f"   - Phonon Hamiltonian: H_ph = sum_{q,j} {w_q} * a_{q,j}^dagger * a_{q,j}")
    H_eph_str = f"g * sum_{q,j} (i*{q_j} / (2*m*{w_q})^(1/2)) * {rho_q} * (a_{q,j} + a_{-q,j}^dagger)"
    print(f"   - Electron-Phonon Interaction: H_e-ph = {H_eph_str}\n")
    
    # --- 3. Integrate out Phonons ---
    print("2. We perform a Gaussian path integral over the phonon fields (a_q, a_q^dagger).")
    print("   This results in an effective interaction term in the action, S_eff, which describes")
    print("   the interaction between electron densities rho_q mediated by a virtual phonon.\n")
    print("   The effective action S_eff is related to the effective interaction V_eff by:")
    print(f"   S_eff = -1/2 * integral(d_tau1 d_tau2) V_eff(q, tau1-tau2) * {rho_q}(tau1) * rho_{-q}(tau2)\n")

    # --- 4. Calculate the Interaction Vertex ---
    # The calculation involves finding the phonon propagator D(q, i*omega_n) and the
    # electron-phonon coupling vertex C_q.
    # The effective interaction V_eff is given by |C_q|^2 * D(q, i*omega_n)
    
    # Coupling vertex constant squared |C_q,j|^2
    C_qj_sq = (g**2 * q_j**2) / (2 * m * w_q)
    
    # Phonon propagator D(q, i*omega_n)
    # This comes from summing the propagators for a and a_dagger:
    # D = -1/(i*omega_n - w_q) - 1/(-i*omega_n - w_q) = 2*w_q / (omega_n^2 + w_q^2)
    D_q = (2 * w_q) / (omega_n**2 + w_q**2)
    
    # The effective interaction for a single polarization 'j'
    V_eff_j = -C_qj_sq * D_q
    
    print("3. For a specific phonon mode (q, j), the effective interaction vertex is V_eff(q, j, i*omega_n):")
    print(f"   V_eff(q, j, i*{omega_n}) = -|C_q,j|^2 * D(q, i*{omega_n})")
    print(f"                       = -({C_qj_sq}) * ({D_q})")
    
    final_V_eff_j = sympy.simplify(V_eff_j)
    print("\n   After simplifying, we get:")
    print(f"   V_eff(q, j, i*{omega_n}) = {final_V_eff_j}\n")
    
    # --- 5. Sum over Polarizations ---
    # The total interaction for a momentum transfer q is the sum over polarizations j.
    # Sum_j (q_j^2) = |q|^2
    q_norm_sq = sympy.Symbol('|q|^2')
    final_V_eff_total = final_V_eff_j.subs(q_j**2, q_norm_sq)

    print("4. To get the total interaction for a momentum transfer q, we sum over all polarizations j:")
    print(f"   V_eff(q, i*{omega_n}) = sum_j V_eff(q, j, i*{omega_n})")
    print("                       Assuming q_j are Cartesian components, sum_j(q_j^2) = |q|^2.")
    print("\n   The final expression for the effective electron-electron interaction is:")
    
    # Using python's f-string formatting to highlight the structure of the final answer
    # Note: sympy.pretty() could be used for fancier printing if needed.
    final_expression = f"-({g**2} * {q_norm_sq}) / ({m} * ({omega_n**2} + {w_q**2}))"
    print(f"\n   V_eff(q, i*{omega_n}) = {final_expression}")

    # --- 6. Static Limit ---
    print("\n   In the static limit (omega_n -> 0), often used in BCS theory, the interaction becomes:")
    static_V_eff = final_V_eff_total.subs(omega_n, 0)
    static_expression = f"-({g**2} * {q_norm_sq}) / ({m} * {w_q**2})"
    print(f"   V_eff(q) = {static_expression}")
    print("\n   This interaction is attractive (note the negative sign), which can lead to the")
    print("   formation of Cooper pairs and superconductivity.")


if __name__ == '__main__':
    derive_effective_interaction()
