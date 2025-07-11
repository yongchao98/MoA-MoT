def derive_effective_interaction():
    """
    Prints the step-by-step derivation of the effective electron-electron
    interaction by integrating out phonon fields.
    """
    print("### Derivation of the Effective Electron-Electron Interaction ###")
    print("\nStep 1: Define the System")
    print("--------------------------")
    print("We start with the total Hamiltonian H = H_el + H_ph + H_el-ph.")
    print("The phonon Hamiltonian is: H_ph = sum_{q,j} w_q * a_{q,j}^dagger * a_{q,j}")
    print("The electron-phonon interaction is: H_el-ph = sum_{q,j} M_{q,j} * n_q * (a_{q,j} + a_{-q,j}^dagger)")
    print("where the coupling M_{q,j} is given by: M_{q,j} = g * i * q_j / (2 * m * w_q)^(1/2)")
    print("and n_q is the electron density operator for momentum q.")

    print("\nStep 2: Path Integral Formulation")
    print("---------------------------------")
    print("The partition function Z is written as a path integral over electron (c) and phonon (a) fields:")
    print("Z = integral D[c*,c,a*,a] exp(-S[c,a])")
    print("Our goal is to integrate out the phonon fields a*, a to find an effective action S_eff[c] for the electrons.")

    print("\nStep 3: Gaussian Integration over Phonon Fields")
    print("------------------------------------------------")
    print("The part of the action involving phonons is quadratic, which allows for an exact Gaussian integral.")
    print("S_ph + S_el-ph = integral d(tau) [ sum_{q,j} a_{q,j}^*(d/d(tau) + w_q)a_{q,j} + H_el-ph ]")
    print("Integrating out the phonon fields (a*, a) yields a new term in the electron action, Delta_S_eff.")
    print("This effective action can be calculated exactly. To second order, it is given by:")
    print("Delta_S_eff = -1/2 * <(S_int)^2>_0")
    print("where <...>_0 denotes an average over the free phonon action and S_int = integral d(tau) H_el-ph.")

    print("\nStep 4: Calculating the Effective Action")
    print("-----------------------------------------")
    print("We compute the average:")
    print("<(S_int)^2>_0 = integral d(tau)d(tau') sum_{q,j,q',j'} M_{q,j}*M_{q',j'}*n_q(tau)*n_{q'}(tau') * <(a_{q,j}(tau)+a_{-q,j}^*(tau))(a_{q',j'}(tau')+a_{-q',j'}^*(tau'))>_0")
    print("The expectation value is non-zero only if q' = -q and j' = j. This gives the phonon propagator D_ph(q, j, tau-tau').")
    print("Delta_S_eff = -1/2 * integral d(tau)d(tau') sum_{q,j} M_{q,j}*M_{-q,j}*n_q(tau)*n_{-q}(tau')*D_ph(q,j,tau-tau')")

    print("\nStep 5: Transforming to Frequency Space")
    print("-----------------------------------------")
    print("In Matsubara frequency space (w_n), the convolution becomes a product:")
    print("Delta_S_eff = -1/(2*beta) * sum_{q,j,n} M_{q,j}*M_{-q,j}*n_{q,n}*n_{-q,-n}*D_ph(q,j,i*w_n)")
    print("We use the following expressions:")
    print("  1. Coupling product: M_{q,j} * M_{-q,j} = (g*i*q_j/...)* (g*i*(-q_j)/...) = g^2 * q_j^2 / (2 * m * w_q)")
    print("  2. Phonon propagator: D_ph(q, j, i*w_n) = 2 * w_q / (w_n^2 + w_q^2)")
    print("  3. Density correlation: n_{-q,-n} = n_{q,n}^*")

    print("\nStep 6: The Effective Interaction Potential")
    print("--------------------------------------------")
    print("Substituting these into the effective action gives:")
    print("Delta_S_eff = -1/(2*beta) * sum_{q,j,n} [g^2 * q_j^2 / (2*m*w_q)] * [2*w_q / (w_n^2 + w_q^2)] * |n_{q,n}|^2")
    print("Delta_S_eff = -1/(2*beta) * sum_{q,j,n} [g^2 * q_j^2 / (m * (w_n^2 + w_q^2))] * |n_{q,n}|^2")
    print("This action corresponds to an effective Hamiltonian H_eff = 1/2 * sum V_eff * n_q^* * n_q.")
    print("The frequency-dependent interaction is: V_eff(q, j, i*w_n) = - g^2 * q_j^2 / (m * (w_n^2 + w_q^2))")

    print("\nStep 7: The Static Interaction")
    print("------------------------------")
    print("The effective electron-electron interaction is typically considered in the static limit, where the energy transfer is zero (w_n -> 0).")
    print("V_eff(q, j) = V_eff(q, j, 0) = - g^2 * q_j^2 / (m * w_q^2)")
    print("This is an attractive interaction (V_eff < 0), which can lead to phenomena like superconductivity.")

    print("\n### Final Result ###")
    print("The effective electron-electron interaction for a given mode (q, j) is:")
    print("V_eff(q, j) = - (g^2 * q_j^2) / (m * w_q^2)")
    print("\nWhere:")
    print("  g: electron-phonon coupling constant")
    print("  q_j: j-th component of the phonon wavevector q")
    print("  m: electron mass")
    print("  w_q: phonon frequency for wavevector q")

if __name__ == '__main__':
    derive_effective_interaction()