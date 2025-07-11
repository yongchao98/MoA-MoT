import sympy

def solve_interaction():
    """
    This function symbolically derives the effective electron-electron interaction
    by integrating out phonon fields, following the provided Hamiltonian.
    """
    # Define the symbolic variables used in the expressions.
    # g: coupling constant, m: mass, w_q: phonon frequency at wavevector q,
    # q_j: j-th component of wavevector q, omega_n: Matsubara frequency for bosons.
    # n_q, n_-q: electron density operators in Fourier space.
    g, m = sympy.symbols('g m', real=True, positive=True)
    w_q = sympy.Symbol('w_q', real=True, positive=True)
    q_j = sympy.Symbol('q_j', real=True)
    omega_n = sympy.Symbol('omega_n', real=True)
    n_q, n_mq = sympy.symbols('n_q n_-q')
    i = sympy.I

    print("Derivation of the Effective Electron-Electron Interaction\n")

    # Step 1: Define the coupling constant M_qj from the interaction Hamiltonian.
    # The interaction is of the form H_int = sum_{q,j} M_{q,j} * n_q * (a_{q,j} + a_{-q,j}^dagger)
    M_qj = g * i * q_j / sympy.sqrt(2 * m * w_q)
    
    # And the corresponding coupling M_{-q,j} for the mode with wavevector -q
    M_mqj = g * i * (-q_j) / sympy.sqrt(2 * m * w_q)

    # Step 2: The strength of the interaction mediated by mode (q,j) is determined
    # by the product M_{q,j} * M_{-q,j}.
    M_prod_j = sympy.simplify(M_qj * M_mqj)

    # Step 3: The free phonon propagator D_0(q, omega_n) in Matsubara space is
    # D_0 = -2*w_q / (omega_n^2 + w_q^2).
    D_0 = -2 * w_q / (omega_n**2 + w_q**2)

    # Step 4: The effective interaction potential V_eff,j from a single mode (q,j) is given by
    # V_eff,j = -M_{q,j}*M_{-q,j}*D_0(q, omega_n).
    V_eff_j = sympy.simplify(-M_prod_j * D_0)

    # Step 5: The total effective potential for a given q is the sum over all components j.
    # This corresponds to replacing sum_j(q_j^2) with the squared magnitude of the vector q, |q|^2.
    q_vec_sq = sympy.Symbol('|q|^2', real=True, positive=True)
    V_eff_total = V_eff_j.subs(q_j**2, q_vec_sq)

    # Step 6: The full effective interaction Hamiltonian term is H_eff = V_eff * n_q * n_-q.
    H_eff = V_eff_total * n_q * n_mq
    
    print("The effective electron-electron interaction Hamiltonian for a given wavevector q is:")
    print(f"H_eff(q, \u03C9_n) = V_eff(q, \u03C9_n) * n_q * n_-q")
    print("\nWhere the effective potential V_eff is:\n")
    
    # Output the components of the final equation as requested.
    # Using specific symbols for clarity in the final print.
    g_s, m_s, w_q_s, q_s, wn_s = sympy.symbols('g m w_q |q| omega_n')
    
    term1_num = g_s**2 * q_s**2
    term1_den = m_s
    term1 = term1_num / term1_den

    term2_num = 1
    term2_den = wn_s**2 + w_q_s**2
    term2 = term2_num / term2_den
    
    final_V_expr = sympy.simplify(V_eff_total.subs({g:g_s, m:m_s, w_q: w_q_s, q_vec_sq: q_s**2, omega_n: wn_s}))
    
    print(f"V_eff(q, \u03C9_n) = ({term1}) * ({term2})")
    
    print("\nThis can be written as the final equation:")
    # Pretty print the final symbolic equation
    final_eq = sympy.Eq(sympy.Symbol("V_eff(q, \u03C9_n)"), final_V_expr)
    sympy.pprint(final_eq, use_unicode=True)
    
    print("\nThus, the full interaction term is:")
    print(f"H_eff(q, \u03C9_n) = {final_V_expr} * n_q * n_-q")

if __name__ == '__main__':
    solve_interaction()