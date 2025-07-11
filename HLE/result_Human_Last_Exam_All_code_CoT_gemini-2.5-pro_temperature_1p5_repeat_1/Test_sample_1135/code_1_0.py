def solve_ccsd_projection_question():
    """
    Analyzes the structure of the CCSD equations to determine for which excited
    Slater determinants the matrix element <Phi_I | H_bar | Phi> is identically zero.
    """
    print("### Analysis of the CCSD Matrix Element <Φ_I | H_bar | Φ> ###\n")

    # 1. Define the fundamental properties of the operators
    H_body_count = 2
    # In CCSD, the cluster operator T = T1 + T2. T2 is the highest excitation operator.
    T_max_excitation_level = 2

    print(f"The electronic Hamiltonian (H) is a {H_body_count}-body operator.")
    print(f"In CCSD, the cluster operator (T) contains up to {T_max_excitation_level}-electron excitations (T1 and T2).\n")

    # 2. Explain the structure of the key state vector H_bar | Phi >
    print("The matrix element <Φ_I | H_bar | Φ> is non-zero only if H_bar|Φ> contains a component of |Φ_I>.")
    print("In Coupled Cluster theory, H_bar |Φ> = (H * exp(T))_connected |Φ>.\n")
    print("We need to find the highest excitation level created when (H * exp(T))_c acts on the reference |Φ>.")

    # 3. The Connectivity Argument
    # A k-body operator can link at most 2k fermion lines.
    # The crucial insight is that a 2-body H cannot connect three or more separate cluster operators (like T1, T2).
    # This limits the excitation level that can be generated.
    # The highest connected term in H_bar |Phi> comes from H acting on a state created by the T operators.
    # For example, H acting on a double excitation (from T2) can create at most a quadruple excitation.
    # H * T2 |Φ> --> generates up to quadruples.
    # H * T2 * T2 |Φ> requires H to connect the two T2 operators. The result is also at most a quadruple.
    max_excitation_in_Hbar_Phi = 4

    print("Due to the two-body nature of the Hamiltonian H, it can be shown that it cannot connect")
    print("enough cluster operators (T1, T2) simultaneously to generate excitations beyond a certain level.")
    print(f"The highest possible excitation level in the expansion of H_bar|Φ> for CCSD is {max_excitation_in_Hbar_Phi} (Quadruples).\n")

    # 4. Conclusion
    print("This means the expansion of H_bar|Φ> looks like:")
    print("H_bar|Φ> = c_0|Φ> + c_S|Φ_S> + c_D|Φ_D> + c_T|Φ_T> + c_Q|Φ_Q>\n")

    print("By the rules of CCSD, the coefficients for singles and doubles are set to zero to solve for the amplitudes:")
    print("<Φ_S | H_bar | Φ> = 0  (by definition)")
    print("<Φ_D | H_bar | Φ> = 0  (by definition)\n")

    print("The projections onto Triples and Quadruples are generally NON-ZERO.")
    print("<Φ_T | H_bar | Φ> ≠ 0")
    print("<Φ_Q | H_bar | Φ> ≠ 0\n")

    print("Because the expansion of H_bar|Φ> stops at quadruples, any projection onto a higher-level determinant")
    print("is identically zero by orthogonality.\n")

    identically_zero_level = max_excitation_in_Hbar_Phi + 1
    print("Therefore, the 'other' excited Slater determinants for which the matrix element is identically zero are those with an excitation level of {} or higher.".format(identically_zero_level))

    final_answer = "Quintuply-excited determinants and all higher-level excited determinants"
    print("\n---")
    print("Final Answer:")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    solve_ccsd_projection_question()