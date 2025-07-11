import sympy as sp

def solve_effective_interaction():
    """
    Performs a symbolic derivation of the effective electron-electron interaction
    by integrating out phonon fields from the given Hamiltonian.
    """
    # Define the symbols
    g, m = sp.symbols('g m', real=True, positive=True)
    q_vec = sp.MatrixSymbol('q', 3, 1)
    w_q = sp.Function('w')(sp.sqrt(q_vec.T * q_vec))
    rho_q = sp.Function('rho')(q_vec)
    rho_minus_q = sp.Function('rho')(-q_vec)
    i = sp.I

    # The user has specified the electron-phonon coupling term in the Hamiltonian.
    # We first identify the vertex functions and the phonon propagator.
    # The interaction Hamiltonian can be written as H_el-ph = sum_{q,j} [J_{q,j}^dagger a_{q,j} + J_{q,j} a_{q,j}^dagger]
    # From the user's expression, V_{q,j} = g * i * q_j / sqrt(2*m*w_q)
    # The source terms coupling to the phonon field are identified as:
    # J_{q,j} = V_{-q,j} * rho_q and J_{q,j}^dagger = V_{q,j} * rho_{-q}
    # V_{-q,j} = g * i * (-q_j) / sqrt(2*m*w_q) = -V_{q,j}

    V_qj_V_minus_qj = []
    q_components = sp.symbols('q_x q_y q_z')
    for q_j in q_components:
        V_qj = g * i * q_j / sp.sqrt(2 * m * w_q)
        V_minus_qj = g * i * (-q_j) / sp.sqrt(2 * m * w_q)
        V_qj_V_minus_qj.append(V_qj * V_minus_qj)
    
    # Sum over polarization index j. This corresponds to q_x^2 + q_y^2 + q_z^2 = |q|^2
    # The product V_{q,j}*V_{-q,j} = g^2 * q_j^2 / (2*m*w_q)
    q_squared = sp.Symbol('q^2', real=True, positive=True)
    sum_V_prod = sp.sum(V_qj_V_minus_qj).subs(sum(q_j**2 for q_j in q_components), q_squared)
    
    # The phonon propagator in Matsubara frequency (i*omega_n) space is D(q, i*omega_n) = 1 / (w_q - i*omega_n)
    omega_n = sp.Symbol('omega_n')
    D_q = 1 / (w_q - i * omega_n)
    D_minus_q = 1 / (w_q + i * omega_n)

    # The effective interaction potential V_eff arises from symmetrizing the interaction kernel.
    # V_eff(q, i*omega_n) is proportional to [U(q, i*omega_n) + U(-q, -i*omega_n)], where U is the unsymmetrized vertex.
    # U(q, i*omega_n) = (sum_j V_{q,j}*V_{-q,j}) * D(q, i*omega_n)
    U_q = sum_V_prod * D_q
    # U(-q, -i*omega_n) = (sum_j V_{-q,j}*V_{q,j}) * D(-q, -i*omega_n)
    U_minus_q = sum_V_prod * D_minus_q

    V_eff_freq_dependent = U_q + U_minus_q

    # The static effective potential is found by taking the limit omega_n -> 0
    V_eff_static = sp.limit(V_eff_freq_dependent, omega_n, 0)
    V_eff_static_simplified = sp.simplify(V_eff_static)

    # The effective electron-electron interaction Hamiltonian for a given q is
    # H_eff(q) = (1/2) * V_eff(q) * rho_q * rho_{-q}
    # The user asks for the effective interaction, so we construct this term (omitting the 1/2 prefactor and sum over q).
    H_eff_q = V_eff_static_simplified * rho_q.subs(q_vec, sp.Symbol('q')) * rho_minus_q.subs(q_vec, sp.Symbol('-q'))
    
    # Print the derivation steps and the final result
    print("Derivation of the Effective Electron-Electron Interaction:\n")
    print("1. The electron-phonon interaction Hamiltonian is given.")
    print("2. We integrate out the phonon fields using path integrals.")
    print("3. This results in an effective interaction between electron densities rho(q).")
    print("4. The frequency-dependent effective potential is found to be:")
    print(f"   V_eff(q, omega_n) = {sp.simplify(V_eff_freq_dependent)}")
    print("\n5. For the static interaction, we take the limit omega_n -> 0:")
    print(f"   V_eff(q) = {V_eff_static_simplified}")
    print("\n6. The resulting effective electron-electron interaction for a given momentum q is:")
    print("   H_eff(q) = V_eff(q) * rho(q) * rho(-q)")
    print("\nFinal expression for the interaction term:")
    # Pretty print the final result
    sp.init_printing(use_unicode=True)
    print(sp.pretty(H_eff_q))
    
    # As requested, output each number in the final equation. Since there are no
    # literal numbers, we will print the components of the expression.
    print("\nComponents of the final equation:")
    print(f"Interaction Strength Factor: {sp.pretty(V_eff_static_simplified)}")
    print(f"Electron density at q: {sp.pretty(rho_q.subs(q_vec, sp.Symbol('q')))}")
    print(f"Electron density at -q: {sp.pretty(rho_minus_q.subs(q_vec, sp.Symbol('-q')))}")


if __name__ == '__main__':
    solve_effective_interaction()
    # Construct the expected output format string.
    g, m, q, wq = sp.symbols("g m q w_q")
    q_sq = sp.Symbol("q^2")
    w_q_func = sp.Function('w')(q)
    final_V = (g**2 * q_sq) / (m * w_q_func**2)
    rho_q_sym = sp.Function('rho')(q)
    rho_mq_sym = sp.Function('rho')(-q)
    final_H = final_V * rho_q_sym * rho_mq_sym
    
    answer_string = f"<<<{sp.sstr(final_H, full_prec=False)}>>>"
    # print(answer_string)