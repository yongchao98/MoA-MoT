import sympy as sp

def derive_effective_interaction():
    """
    This function symbolically derives the effective electron-electron interaction
    mediated by phonons by integrating out the phonon degrees of freedom.
    """
    # 1. Define the symbolic variables used in the problem.
    # g: coupling constant, m: mass
    g, m = sp.symbols('g m', real=True, positive=True)
    # q_j: j-th component of wavevector, w_q: phonon frequency, nu_n: Matsubara frequency
    q_j, w_q, nu_n = sp.symbols('q_j w_q nu_n', real=True)
    # n_q, n_{-q}: electron density operators in momentum space
    n_q = sp.Symbol('n_q')
    n_minus_q = sp.Symbol('n_{-q}')
    # |q|^2: squared magnitude of the wavevector q
    q_vec_sq = sp.Symbol('|q|^2')

    # 2. Define the electron-phonon coupling constant from the Hamiltonian.
    # The interaction term is H_int = sum_{q,j} C_{q,j} * n_q * (a_{q,j} + a_{-q,j}^dagger)
    # with the vertex factor C_{q,j} given by:
    C_qj = g * sp.I * q_j / sp.sqrt(2 * m * w_q)

    # 3. Define the phonon propagator in the imaginary frequency domain.
    # The interaction involves the phonon displacement field, which is proportional to (a_q + a_{-q}^dagger).
    # The propagator for this field in imaginary time formalism is D(q, i*nu_n) = 2*w_q / (w_q^2 + nu_n^2).
    D_q_nun = 2 * w_q / (w_q**2 + nu_n**2)

    # 4. Calculate the effective interaction potential V_eff.
    # This is obtained from second-order perturbation theory, which is equivalent to the
    # path integral result. The potential is V_eff = -C_{q,j} * C_{-q,j} * D(q, i*nu_n).
    # C_{-q,j} is found by replacing q_j with -q_j.
    C_minus_qj = C_qj.subs(q_j, -q_j)
    
    # Calculate the potential for a single polarization j.
    V_eff_j = -C_qj * C_minus_qj * D_q_nun
    V_eff_j_simplified = sp.simplify(V_eff_j)

    # The total potential is the sum over all polarizations j.
    # Assuming an isotropic system, we can replace the sum sum_j(q_j^2) with |q|^2.
    V_eff_total = V_eff_j_simplified.subs(q_j**2, q_vec_sq)

    # 5. Construct the effective interaction term in the action for a single q.
    # The effective interaction term is Delta_S_q = (1/2) * V_eff * n_q * n_{-q}.
    # (The summations over frequency and polarization are part of V_eff here).
    H_eff_q = sp.Rational(1, 2) * V_eff_total * n_q * n_minus_q

    # 6. Print the results clearly.
    print("The effective electron-electron interaction is derived by integrating out the phonon fields.")
    print("The process yields an effective potential that describes the interaction between electron density fluctuations.")
    print("-" * 70)
    print("Key components of the derivation:")
    print(f"1. Electron-Phonon Coupling Vertex C(q,j): {C_qj}")
    print(f"2. Phonon Propagator D(q, i*nu_n): {D_q_nun}")
    print(f"3. Effective Potential V_eff(q, i*nu_n) after summing over polarizations j:")
    # Using pretty print for better formatting of the expression
    sp.pprint(V_eff_total, use_unicode=False)
    print("-" * 70)

    print("\nThe final effective electron-electron interaction term for a given wavevector q is:")
    print("H_eff(q) =")
    sp.pprint(H_eff_q, use_unicode=False)
    print("\nThis expression represents an attractive interaction (due to the minus sign) between electron density")
    print("fluctuations, mediated by the exchange of virtual phonons.")
    print("\nEach symbol in the final equation represents:")
    print(f"  {g}: electron-phonon coupling strength")
    print(f"  {m}: electron mass")
    print(f"  {q_vec_sq}: the squared magnitude of the wavevector q")
    print(f"  {w_q}: the phonon frequency at wavevector q")
    print(f"  {nu_n}: the bosonic Matsubara frequency, indicating the energy transferred in the interaction")
    print(f"  {n_q}, {n_minus_q}: the electron density operators for wavevectors q and -q")


if __name__ == '__main__':
    derive_effective_interaction()
