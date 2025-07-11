# Define the truth values. We'll use strings for clarity.
T, G, F = "T", "G", "F"
vals = [T, G, F]

# Define the truth tables for the 3-valued logic LP (Logic of Paradox)
neg_table = {T: F, G: G, F: T}
conj_table = {
    T: {T: T, G: G, F: F},
    G: {T: G, G: G, F: F},
    F: {T: F, G: F, F: F},
}
disj_table = {
    T: {T: T, G: T, F: T},
    G: {T: T, G: G, F: G},
    F: {T: T, G: G, F: F},
}

def imply(a, b):
    """Computes the value of a -> b as not(a) or b."""
    return disj_table[neg_table[a]][b]

def check_argument_K():
    """
    Checks the validity of argument K: A & B |- (¬A ∨ ¬B) -> (A & B)
    using t-entailment (preservation of T).
    """
    print("Checking validity of argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
    print("Using t-entailment: if premise is T, conclusion must be T.\n")

    is_valid = True
    
    # We only need to check the cases where the premise is T.
    # The premise A ∧ B is T only when A=T and B=T.
    v_A = T
    v_B = T

    print(f"Testing the only case where the premise is T: v(A)={v_A}, v(B)={v_B}")

    # Evaluate premise
    premise_val = conj_table[v_A][v_B]
    print(f"  Premise: v(A ∧ B) = v({v_A} ∧ {v_B}) = {premise_val}")
    
    # If premise is T, check conclusion
    if premise_val == T:
        # Evaluate conclusion: (¬A ∨ ¬B) → (A ∧ B)
        neg_A = neg_table[v_A]
        neg_B = neg_table[v_B]
        print(f"    Evaluating conclusion part 1: v(¬A ∨ ¬B)")
        print(f"      v(¬{v_A}) = {neg_A}")
        print(f"      v(¬{v_B}) = {neg_B}")
        
        part1_val = disj_table[neg_A][neg_B]
        print(f"      v(¬A ∨ ¬B) = v({neg_A} ∨ {neg_B}) = {part1_val}")

        print(f"    Evaluating conclusion part 2: v(A ∧ B)")
        part2_val = conj_table[v_A][v_B]
        print(f"      v(A ∧ B) = v({v_A} ∧ {v_B}) = {part2_val}")

        conclusion_val = imply(part1_val, part2_val)
        print(f"    Final Conclusion Value: v(({part1_val}) → ({part2_val})) = {conclusion_val}\n")

        if conclusion_val != T:
            is_valid = False
            print(f"  COUNTEREXAMPLE FOUND: Premise is T, but conclusion is {conclusion_val}.")
    
    if is_valid:
        print("The argument is valid under t-entailment.")
    else:
        print("The argument is invalid under t-entailment.")

check_argument_K()