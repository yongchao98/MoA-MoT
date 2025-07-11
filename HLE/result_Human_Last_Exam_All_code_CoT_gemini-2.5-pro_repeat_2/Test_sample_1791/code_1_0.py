import sympy as sp

def derive_effective_interaction():
    """
    This function symbolically derives the effective electron-electron interaction
    by integrating out the phonon fields from the given Hamiltonian.
    It prints the steps of the derivation and the final equation.
    """
    # Define the symbolic variables used in the expressions.
    # g: electron-phonon coupling strength
    # m: mass parameter from the coupling
    # q_j: j-th component of the phonon momentum vector q
    # w_q: phonon frequency for momentum q
    # n_q, n_mq: electron density operators for momentum q and -q
    g, m = sp.symbols('g m', real=True, positive=True)
    q_j = sp.Symbol('q_j', real=True)
    w_q = sp.Symbol('w_q', real=True, positive=True)
    n_q, n_mq = sp.symbols('n_q n_-q')
    H_eff_int = sp.Symbol('H_eff_int(q,j)')

    print("Derivation of the Effective Electron-Electron Interaction\n")

    print("Step 1: Identify the fundamental electron-phonon vertex V_qj.")
    # The interaction Hamiltonian is H_int = g * sum(...) * (i * q_j / sqrt(2*m*w_q)) * n_q * (a_{q,j} + a_{-q,j}^dagger).
    # In the path integral formalism, this interaction acts as a source for the phonon field.
    # The fundamental vertex V_qj connects an electron density fluctuation to a phonon.
    V_qj_expr = g * sp.I * q_j / sp.sqrt(2 * m * w_q)
    print(f"The vertex is defined from the coupling as: V_qj = {V_qj_expr}\n")

    print("Step 2: Calculate the squared magnitude of the vertex, |V_qj|^2.")
    # The strength of the resulting interaction is proportional to |V_qj|^2.
    V_qj_mag_sq = sp.simplify(sp.Abs(V_qj_expr)**2)
    print(f"|V_qj|^2 = g^2 * q_j^2 / (2 * m * w_q)")
    print(f"Simplified: |V_qj|^2 = {V_qj_mag_sq}\n")

    print("Step 3: Determine the static effective electron-electron potential V_eff(q,j).")
    # Integrating out the phonon fields results in an effective potential. In the static limit (frequency omega -> 0),
    # the standard result from the path integral is V_eff(q,j) = -2 * |V_qj|^2 / w_q.
    V_eff_expr = -2 * V_qj_mag_sq / w_q
    print(f"The potential is given by V_eff(q,j) = -2 * |V_qj|^2 / w_q.")
    print(f"Substituting the expression for |V_qj|^2:")
    print(f"V_eff(q,j) = -2 * ({V_qj_mag_sq}) / {w_q} = {V_eff_expr}")
    print("The negative sign indicates that the phonon-mediated interaction is attractive.\n")

    print("Step 4: Construct the effective interaction Hamiltonian term.")
    # The interaction is of the form H_eff_int = (1/2) * sum_q V_eff(q) * n_q * n_-q.
    # The factor of 1/2 is a convention to avoid double-counting when summing over all momenta q.
    # As requested, we provide the term for a single mode (q,j) before summation.
    final_hamiltonian_expr = (sp.S(1)/2) * V_eff_expr * n_q * n_mq

    print("The effective interaction term for mode (q,j) is H_eff_int(q,j) = (1/2) * V_eff(q,j) * n_q * n_-q.")
    print("\n--- FINAL EQUATION ---")
    print("The final equation for the effective electron-electron interaction, showing each component, is:")

    # Explicitly print each "number" or symbol in the final equation.
    # The coefficient is -g**2*q_j**2/(2*m*w_q**2).
    print(f"{H_eff_int} = - ( 1 / 2 ) * ( g**2 * q_j**2 ) / ( m * w_q**2 ) * {n_q} * {n_mq}")

derive_effective_interaction()