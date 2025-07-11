def derive_effective_interaction():
    """
    This script symbolically derives the effective electron-electron interaction
    mediated by phonons by performing a path integral over the phonon fields.
    """

    print("### Derivation of the Effective Electron-Electron Interaction ###\n")

    # Step 1: Define the Hamiltonians
    print("Step 1: Define the relevant parts of the Hamiltonian.")
    print("------------------------------------------------------")
    print("The phonon Hamiltonian is given by:")
    print("H_ph = sum_{q, j} w_q * a_{q,j}^dagger * a_{q,j}")
    print("\nThe electron-phonon interaction Hamiltonian is:")
    print("H_e-ph = sum_{q, j} V_{q,j} * rho_q * (a_{q,j} + a_{-q,j}^dagger)")
    print("where rho_q is the electron density operator and the coupling V_{q,j} is:")
    print("V_{q,j} = g * i * q_j / (2 * m * w_q)^(1/2)\n")

    # Step 2: Formulate the effective action
    print("Step 2: Formulate the effective action.")
    print("---------------------------------------")
    print("Integrating out the phonon fields in the path integral generates an effective action for the electrons.")
    print("To second order in the coupling, this effective action term, Delta_S, is given by:")
    print("Delta_S = (1/2) * integral[d(tau) d(tau')] <T_tau H_e-ph(tau) * H_e-ph(tau')>_0")
    print("where the expectation value is taken with respect to the free phonon Hamiltonian.\n")
    print("This corresponds to an effective interaction term in the action of the form:")
    print("Delta_S = (1/2) * sum_{q, n} V_eff(q, i*nu_n) * rho_{q,-n} * rho_{-q,n}\n")

    # Step 3: Evaluate the expectation value
    print("Step 3: Evaluate the expectation value using the phonon propagator.")
    print("-----------------------------------------------------------------")
    print("The expectation value only connects terms with opposite momenta (q' = -q).")
    print("This leads to:")
    print("Delta_S = (1/2) * sum_{q,j,n} V_{q,j} * V_{-q,j} * D(q,j; i*nu_n) * rho_{q,-n} * rho_{-q,n}")
    print("where D(q,j; i*nu_n) is the Fourier transform of the free phonon propagator:")
    print("D(q,j; tau-tau') = <T_tau (a_{q,j} + a_{-q,j}^dagger)(tau) * (a_{-q,j} + a_{q,j}^dagger)(tau')>_0")
    print("\nThe value of the propagator in Matsubara frequency (nu_n) space is:")
    propagator_expr = "-2 * w_q / (w_q^2 + nu_n^2)"
    print(f"D(q,j; i*nu_n) = {propagator_expr}\n")

    # Step 4: Substitute the coupling constants
    print("Step 4: Substitute the specific form of the coupling V_{q,j}.")
    print("-------------------------------------------------------------")
    print("First, we compute the product V_{q,j} * V_{-q,j}:")
    print("V_{q,j} = g * i * q_j / sqrt(2*m*w_q)")
    print("V_{-q,j} = g * i * (-q_j) / sqrt(2*m*w_q) = -V_{q,j}")
    print("So, V_{q,j} * V_{-q,j} = - (V_{q,j})^2 = - (g * i * q_j / sqrt(2*m*w_q))^2")
    v_prod_expr = "g^2 * q_j^2 / (2 * m * w_q)"
    print(f"V_{q,j} * V_{-q,j} = {v_prod_expr}\n")

    # Step 5: Combine terms to find the effective potential for a given q
    print("Step 5: Combine terms to find the effective potential V_eff.")
    print("-----------------------------------------------------------")
    print("The effective potential V_eff(q, i*nu_n) is the coefficient of (1/2) * rho_{q,-n} * rho_{-q,n}.")
    print("V_eff(q, i*nu_n) = sum_j V_{q,j} * V_{-q,j} * D(q,j; i*nu_n)")
    print("Substituting the expressions from the previous steps:")
    print(f"V_eff(q, i*nu_n) = sum_j [ {v_prod_expr} ] * [ {propagator_expr} ]")
    print("\nSimplifying the expression by cancelling terms:")
    intermediate_expr = "- sum_j g^2 * q_j^2 / m"
    print(f"V_eff(q, i*nu_n) = ( {intermediate_expr} ) * ( 1 / (w_q^2 + nu_n^2) )")
    print("\nFinally, performing the sum over polarization index j (e.g., j=x,y,z):")
    print("sum_j q_j^2 = q_x^2 + q_y^2 + q_z^2 = |q|^2 = q^2")
    
    # Final Result
    print("\n### Final Result ###")
    print("The effective electron-electron interaction potential for a given momentum q is:")
    final_potential_term1 = "-g^2"
    final_potential_term2 = "q^2"
    final_potential_term3 = "m"
    final_potential_term4 = "1 / (w_q^2 + nu_n^2)"
    print("V_eff(q, i*nu_n) = ( {} * {} / {} ) * ( {} )".format(
        final_potential_term1,
        final_potential_term2,
        final_potential_term3,
        final_potential_term4
    ))
    print("\nThis represents an attractive interaction between electrons mediated by phonons.")


if __name__ == '__main__':
    derive_effective_interaction()
    
# The final result is the effective potential V_eff(q, i*nu_n).
# Let's extract the components as requested for the final answer format.
g_term = "-g^2"
q_term = "q^2"
m_term = "m"
freq_term = "1/(w_q^2 + nu_n^2)"

# The user wants the answer in the format <<<answer content>>>.
# The most reasonable "answer content" is the final expression for V_eff.
# The instruction "output each number in the final equation" is interpreted as showing the components,
# which the script does. The final expression itself is the complete answer.
final_answer = "V_eff(q, i*nu_n) = -g^2 * q^2 / m * (1 / (w_q^2 + nu_n^2))"
# The prompt is a bit ambiguous about whether the code output IS the answer, or if I should put
# the final string into the <<<>>> format. I'll stick to the prompt's instruction to provide code.
# The code itself, when run, will print the final answer in a clear way.
