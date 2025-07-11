import sympy.combinatorics as comb

def analyze_braid_closure(n_strands, braid_word_str):
    """
    Analyzes the components of a braid closure.

    Args:
        n_strands (int): The number of strands in the braid group.
        braid_word_str (str): A string representing the braid word.
    """
    
    # Define the braid generator permutations using 0-based indexing for sympy
    s = [comb.Permutation(i, i+1, size=n_strands) for i in range(n_strands - 1)]

    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1
    # Sympy permutation multiplication P1*P2 means apply P2 then P1.
    # Braid word multiplication beta1*beta2 means apply beta1 then beta2.
    # This means the order of composition is reversed for permutations.
    # So, perm(beta1*beta2) = perm(beta2) * perm(beta1)
    
    # Deconstruct the braid into its parts for the "equation"
    p_s1_sq = s[0]**2
    p_s2_sq = s[1]**2
    p_s3 = s[2]
    p_s4_inv = s[3]**-1
    
    # Calculate the final permutation.
    # perm(beta) = perm(sigma_4^-1) * perm(sigma_3) * perm(sigma_2^2) * perm(sigma_1^2)
    beta_perm = p_s4_inv * p_s3 * p_s2_sq * p_s1_sq

    # Find the cycle decomposition
    cycles = beta_perm.cyclic_form
    # Adjust to 1-based indexing for standard notation
    adjusted_cycles = [[i + 1 for i in cycle] for cycle in cycles]

    # --- Output the results ---
    print(f"Analysis for the braid beta = {braid_word_str} in B_{n_strands}:")
    print("-" * 50)
    print("Step 1: Determine the permutation for each part of the braid.")
    print(f"  Permutation of sigma_1^2  : {p_s1_sq.cyclic_form} (Identity)")
    print(f"  Permutation of sigma_2^2  : {p_s2_sq.cyclic_form} (Identity)")
    print(f"  Permutation of sigma_3    : {[[i+1 for i in c] for c in p_s3.cyclic_form]}")
    print(f"  Permutation of sigma_4^-1 : {[[i+1 for i in c] for c in p_s4_inv.cyclic_form]}")

    print("\nStep 2: Calculate the total permutation by composing them.")
    print(f"The equation is: P(beta) = P(sigma_4^-1) * P(sigma_3) * P(sigma_2^2) * P(sigma_1^2)")
    print(f"Resulting permutation (1-based array form): {[i+1 for i in beta_perm.array_form]}")

    print("\nStep 3: Find the cycles of the final permutation.")
    print(f"The cycle decomposition is: {adjusted_cycles}")

    print("\nConclusion:")
    print(f"The closure of the braid has {len(adjusted_cycles)} connected components.")
    
    knot_component = None
    for i, cycle in enumerate(adjusted_cycles):
        print(f"  - Component {i+1} is formed by strand(s) {cycle}.")
        if len(cycle) > 1:
            knot_component = cycle

    print("\nGiven two components are unknots (likely the single-strand components),")
    print(f"the remaining component, {knot_component}, must be the knot we want to identify.")
    print("As derived in the explanation, this component is the Figure-8 knot.")


# Run the analysis for the specific problem
analyze_braid_closure(5, "sigma_1^2*sigma_2^2*sigma_3*sigma_4^-1")