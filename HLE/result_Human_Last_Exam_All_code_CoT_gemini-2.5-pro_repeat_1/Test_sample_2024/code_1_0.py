import sys

# Define the set of truth values
TRUTH_VALUES = {0, 0.5, 1} # False, Indeterminate, True

# Define the worlds. In S5, they form an equivalence class.
WORLDS = {'w1', 'w2', 'w3'}

def implication_k3(p, q):
    """Defines the Kleene K3 implication."""
    if p == 0:
        return 1
    elif p == 1:
        return q
    elif p == 0.5:
        if q == 1:
            return 1
        else:
            return 0.5
    raise ValueError("Invalid truth value")

def check_bivalence_of_T():
    """
    Step 1: Prove by contradiction that T(x,y,z) cannot have truth value 0.5.
    The proof assumes the 'Axiom Truth Value' must evaluate to 1.
    Axiom: T(x,y,z) -> Box(ForAll w (R(z,w) -> T(x,y,w)))
    """
    print("Step 1: Analyzing the 'Axiom Truth Value' to determine properties of predicate T.")
    print("Axiom: T(x, y, z) -> □(∀w (R(z, w) -> T(x, y, w)))")
    print("This axiom must have a truth value of 1 in all worlds for any x, y, z.")
    
    # Assume for contradiction that v_wc(T(x,y,z)) = 0.5 for some world wc.
    val_T_in_wc = 0.5
    print(f"\nAssume for contradiction, there is a T-statement, let's call it A, with a truth value of {val_T_in_wc}.")
    
    # Let C be the consequent of the axiom: C = □(...)
    # The axiom's value is implication_k3(v(A), v(C)). This must be 1.
    # For implication_k3(0.5, v(C)) to be 1, v(C) must be 1.
    val_C_in_wc = 1
    print(f"For the axiom '{val_T_in_wc} -> C' to be 1, the consequent C must have a truth value of {val_C_in_wc}.")
    
    # v_wc(C) = v_wc(□(...)) = min(v_w'(...) for w' in WORLDS).
    # If the minimum is 1, the value in every world must be 1.
    # Let's analyze the formula inside the Box: P_axiom = ∀w (R(z,w) -> T(x,y,w))
    print(f"The consequent C being {val_C_in_wc} means the formula inside the □, let's call it P_axiom, must have a truth value of 1 in all worlds.")
    
    # Specifically, v_wc(P_axiom) must be 1.
    # This means for any world w_test, v_wc(R(z,w_test) -> T(x,y,w_test)) must be 1.
    # Since R is reflexive, R(z,z) is true. Let's pick w_test = z.
    # The implication becomes: True -> T(x,y,z).
    # Its value in wc is implication_k3(1, v_wc(T(x,y,z))).
    val_implication_in_axiom = implication_k3(1, val_T_in_wc)
    print(f"Analyzing P_axiom in world wc shows that 'True -> A' must be 1. With v(A)={val_T_in_wc}, this implication evaluates to {val_implication_in_axiom}.")

    if val_implication_in_axiom != 1:
        print(f"\nThis is a contradiction! We concluded P_axiom must be 1, but found an instance where it evaluates to {val_implication_in_axiom}.")
        print("Therefore, our initial assumption was false.")
        print("Conclusion: No statement T(x,y,z) can have the truth value 0.5. T is bivalent (its value can only be 0 or 1).\n")
        return True # T is bivalent
    else:
        # This path shouldn't be taken with K3 logic.
        print("\nNo contradiction found. The reasoning depends on the specific 3-valued logic used.")
        return False

def evaluate_target_statement(T_is_bivalent):
    """
    Step 2 & 3: Evaluate the target statement using the bivalence of T.
    Statement: □(∀x∀y∀z (T(x, y, z) -> □(T(x, y, z))))
    """
    if not T_is_bivalent:
        print("Cannot proceed without establishing the bivalence of T.")
        return

    print("Step 2: Evaluating the inner formula P = ∀x∀y∀z (T(x, y, z) -> □(T(x, y, z))).")
    print("Since T is bivalent, we only need to check two cases for the truth value of any T-statement, T_stmt.")

    # Case 1: v(T_stmt) = 0
    antecedent_val = 0
    implication_val_case_1 = implication_k3(antecedent_val, 0.5) # Consequent doesn't matter
    print(f"  - Case A: If v(T_stmt) = {antecedent_val}, the implication '{antecedent_val} -> ...' is {implication_val_case_1}.")

    # Case 2: v(T_stmt) = 1 in a world wc.
    antecedent_val = 1
    # As per Step 1, if T(x,y,z) is true in one world, the axiom forces it to be true in all accessible worlds.
    # So v_w'(T(x,y,z)) = 1 for all w'.
    # Therefore, v_wc(□(T(x,y,z))) = min({1, 1, 1}) = 1.
    consequent_val = 1
    implication_val_case_2 = implication_k3(antecedent_val, consequent_val)
    print(f"  - Case B: If v(T_stmt) = {antecedent_val}, the axiom implies its value is 1 in all accessible worlds.")
    print(f"    Therefore, v(□(T_stmt)) = {consequent_val}.")
    print(f"    The implication '{antecedent_val} -> {consequent_val}' is {implication_val_case_2}.")

    if implication_val_case_1 == 1 and implication_val_case_2 == 1:
        val_P = 1
        print(f"\nIn all cases, the implication is 1. Thus, the universally quantified formula P is always {val_P}.\n")
    else:
        print("\nCould not determine the value of P.")
        return

    print("Step 3: Evaluating the final expression □(P).")
    # v_w1(□P) = min(v_w1(P), v_w2(P), v_w3(P))
    final_value = min({val_P, val_P, val_P})
    print(f"The value of □(P) is the minimum of P's value across all accessible worlds.")
    print(f"v_w1(□P) = min({val_P}, {val_P}, {val_P}) = {final_value}")
    
    return final_value

if __name__ == "__main__":
    T_is_bivalent = check_bivalence_of_T()
    final_truth_value = evaluate_target_statement(T_is_bivalent)
    print("\n" + "="*40)
    print(f"The final truth value of the statement is: {final_truth_value}")
    print("="*40)
    # Final answer format
    sys.stdout.write(f'<<<{final_truth_value}>>>')
