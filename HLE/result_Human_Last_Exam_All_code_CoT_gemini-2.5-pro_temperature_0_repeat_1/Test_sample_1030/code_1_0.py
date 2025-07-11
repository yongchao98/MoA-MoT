import itertools

def neg(v):
    """Implements 3-valued negation."""
    if v == 'T': return 'F'
    if v == 'F': return 'T'
    return 'G'  # v == 'G'

def conj(v1, v2):
    """Implements 3-valued conjunction."""
    order = {'F': 0, 'G': 1, 'T': 2}
    rev_order = {0: 'F', 1: 'G', 2: 'T'}
    return rev_order[min(order[v1], order[v2])]

def disj(v1, v2):
    """Implements 3-valued disjunction."""
    order = {'F': 0, 'G': 1, 'T': 2}
    rev_order = {0: 'F', 1: 'G', 2: 'T'}
    return rev_order[max(order[v1], order[v2])]

def impl(v1, v2):
    """Implements 3-valued implication as not(v1) or v2."""
    return disj(neg(v1), v2)

def evaluate_argument_K():
    """
    Checks the validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    in a 3-valued logic where only 'T' is the designated truth value.
    """
    print("="*50)
    print("Analyzing Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print("="*50)
    print("System: 3-valued logic (T, G, F)")
    print("Designated value for truth: {T}\n")

    is_valid = True
    checked_cases = 0
    
    values = ['T', 'G', 'F']
    
    # Iterate through all possible assignments for A and B
    for a_val in values:
        for b_val in values:
            # Evaluate the premise: A ∧ B
            premise_val = conj(a_val, b_val)
            
            # Validity requires that IF the premise is designated (T),
            # THEN the conclusion must also be designated (T).
            if premise_val == 'T':
                checked_cases += 1
                print(f"Found case with designated premise: A='{a_val}', B='{b_val}'")
                print(f"  Premise evaluation: {a_val} ∧ {b_val} = {premise_val} (Designated)")
                
                # Evaluate the conclusion: (¬A ∨ ¬B) → (A ∧ B)
                neg_a = neg(a_val)
                neg_b = neg(b_val)
                disj_part = disj(neg_a, neg_b)
                conj_part = conj(a_val, b_val)
                conclusion_val = impl(disj_part, conj_part)
                
                print("  Conclusion evaluation:")
                print(f"    (¬{a_val} ∨ ¬{b_val}) → ({a_val} ∧ {b_val})")
                print(f"    = ({neg_a} ∨ {neg_b}) → {conj_part}")
                print(f"    = {disj_part} → {conj_part}")
                print(f"    = {conclusion_val}")

                if conclusion_val == 'T':
                    print("  Conclusion is Designated. The argument holds for this case.\n")
                else:
                    print(f"  Conclusion is NOT Designated. The argument is INVALID.\n")
                    is_valid = False
                    break
        if not is_valid:
            break
            
    if checked_cases == 0:
        # This case is vacuously true, but we found a case.
        print("No cases found where the premise is designated. The argument is vacuously valid.")

    print("-" * 50)
    if is_valid:
        print("Result: The argument K is VALID in the specified system.")
    else:
        print("Result: The argument K is INVALID in the specified system.")
    print("-" * 50)

# Run the analysis
evaluate_argument_K()