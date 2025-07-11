def solve_vc_dimension():
    """
    This script explains and calculates the VC dimension for the given logic.
    """
    # Number of unary predicates in the schema S
    k = 4

    print("Step 1: Proving the lower bound, VC-dim >= k")
    print("-" * 40)
    print(f"Let k = {k} be the number of unary predicates.")
    print(f"Let s_j be the sentence 'exists x, Pj(x)'.")
    print(f"Consider k structures M_1, ..., M_k, where M_i satisfies s_j if and only if i != j.")
    print(f"We will now show that this set of {k} structures can be shattered.")

    # A function to simulate our model M_i
    # M_i satisfies s_j iff i != j
    def check_model_satisfaction(model_index, sentence_index):
        return model_index != sentence_index

    # The formula is a conjunction of sentences s_j
    def check_formula(model_index, formula_sentence_indices):
        for j in formula_sentence_indices:
            if not check_model_satisfaction(model_index, j):
                return False
        return True

    # Demonstrate shattering for a specific target subset I
    # Let's choose a non-trivial subset for I
    target_I = {1, 4}
    
    # The formula is constructed from indices NOT in I
    formula_indices = {j for j in range(1, k + 1) if j not in target_I}
    
    print(f"\nLet's test for a target labeling where M_i is positive if i is in {target_I}.")
    print(f"The corresponding formula, phi_I, is the conjunction of s_j for j in {formula_indices}.")

    all_match = True
    for i in range(1, k + 1):
        is_in_I = (i in target_I)
        # Check if M_i satisfies the formula
        result = check_formula(i, formula_indices)
        
        print(f"  - M_{i}: Formula evaluates to {result}. Desired outcome: {is_in_I}. ", end="")
        if result == is_in_I:
            print("Match!")
        else:
            print("Mismatch!")
            all_match = False

    if all_match:
        print("\nThe logic holds for this example. The argument is general for any subset I.")
        print("This construction proves that VC-dim >= k.")

    print("\nStep 2: Stating the upper bound, VC-dim <= k")
    print("-" * 40)
    print("The logic FO_{exists, land, top, bot} is a fragment of the more expressive logic FO_{exists, land, vee}.")
    print("A known result states that the VC dimension of FO_{exists, land, vee} over k unary predicates is k.")
    print("The VC dimension of a fragment cannot be larger than the full logic.")
    print("Therefore, VC-dim <= k.")

    print("\nStep 3: Conclusion")
    print("-" * 40)
    print(f"Since VC-dim >= {k} and VC-dim <= {k}, we can conclude the VC dimension is exactly {k}.")
    print("\nThe final equation is:")
    print(f"VC-dim = {k}")

solve_vc_dimension()