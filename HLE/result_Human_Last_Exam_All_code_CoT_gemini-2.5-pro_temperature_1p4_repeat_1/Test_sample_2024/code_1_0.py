def kleene_implication(a, b):
    """
    Calculates the truth value of 'a implies b' using Kleene's 3-valued logic.
    Truth values are 0 (False), 0.5 (Indeterminate), 1 (True).
    """
    return max(1 - a, b)

def necessity_operator(values_in_accessible_worlds):
    """
    Calculates the truth value of '□A' as the minimum value of A in all accessible worlds.
    """
    return min(values_in_accessible_worlds)

def evaluate_statement():
    """
    Evaluates the truth value of the statement based on the logical derivation.
    """
    print("Based on the logical derivation, for any proposition p = T(x,y,z),")
    print("its truth value must be constant across all accessible worlds and must be either 0 or 1.\n")

    # Case 1: The proposition p = T(x,y,z) is false in the equivalence class.
    val_p_case_1 = 0
    # The truth value of p is the same in w1, w2, w3.
    v_p_in_worlds_1 = [val_p_case_1, val_p_case_1, val_p_case_1]
    # The value of □p is the minimum of these values.
    v_box_p_1 = necessity_operator(v_p_in_worlds_1)
    # The value of the implication p -> □p
    result_1 = kleene_implication(val_p_case_1, v_box_p_1)
    
    print(f"Case 1: v(T(x,y,z)) = {val_p_case_1}")
    print(f"  - The statement T(x,y,z) -> □(T(x,y,z)) becomes {val_p_case_1} -> {v_box_p_1}")
    print(f"  - The truth value is max(1 - {val_p_case_1}, {v_box_p_1}) = {result_1}\n")

    # Case 2: The proposition p = T(x,y,z) is true in the equivalence class.
    val_p_case_2 = 1
    # The truth value of p is the same in w1, w2, w3.
    v_p_in_worlds_2 = [val_p_case_2, val_p_case_2, val_p_case_2]
    # The value of □p is the minimum of these values.
    v_box_p_2 = necessity_operator(v_p_in_worlds_2)
    # The value of the implication p -> □p
    result_2 = kleene_implication(val_p_case_2, v_box_p_2)

    print(f"Case 2: v(T(x,y,z)) = {val_p_case_2}")
    print(f"  - The statement T(x,y,z) -> □(T(x,y,z)) becomes {val_p_case_2} -> {v_box_p_2}")
    print(f"  - The truth value is max(1 - {val_p_case_2}, {v_box_p_2}) = {result_2}\n")

    print("Since the inner formula always evaluates to 1, the universally quantified formula")
    print("F = ∀x∀y∀z(...) also has a value of 1. Consequently, □F must also be 1.")
    
    final_truth_value = min(result_1, result_2) # This is a simplification. The real minimum is over all x,y,z.
    print(f"\nThe determined truth value of the full statement in world w1 is: {final_truth_value}")
    
    return final_truth_value

# Run the evaluation
evaluate_statement()