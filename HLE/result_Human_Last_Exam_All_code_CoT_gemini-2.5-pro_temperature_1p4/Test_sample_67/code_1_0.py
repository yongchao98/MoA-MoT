def check_formulas():
    """
    This script models a Kripke structure with decreasing domains to test
    the Barcan Formula and its Converse.
    """
    # Define the Kripke Model
    # W = Worlds, R = Accessibility Relation, D = Domains, I = Interpretation
    W = {'w', 'w1', 'w2'}
    R = {'w': {'w1', 'w2'}, 'w1': set(), 'w2': set()}
    D = {'w': {'a', 'b'}, 'w1': {'a'}, 'w2': {'b'}}

    # Interpretation for predicate 'phi'
    # I[world] = set of individuals with property phi in that world
    I = {
        'w1': {'a'},
        'w2': {'b'}
    }

    print("--- Model Definition ---")
    print(f"Worlds: {W}")
    print(f"Domains: {D}")
    print("Accessibility: w -> w1, w -> w2")
    print("This is a decreasing domain model because D(w1) and D(w2) are subsets of D(w).")
    print(f"Interpretation of phi: phi(a) is true in w1, phi(b) is true in w2.")
    print("-" * 25)

    # --- Evaluate Barcan Formula: □∃x φ(x) → ∃x □φ(x) ---
    print("\nEvaluating Barcan Formula: □∃x φ(x) → ∃x □φ(x)\n")

    # Antecedent: □∃x φ(x) at w
    # This means for all w' accessible from w, ∃x φ(x) is true in w'.
    antecedent_bf_holds_at_w1 = any(x in I.get('w1', set()) for x in D['w1'])
    antecedent_bf_holds_at_w2 = any(x in I.get('w2', set()) for x in D['w2'])
    antecedent_bf_is_true = antecedent_bf_holds_at_w1 and antecedent_bf_holds_at_w2
    print(f"1. Antecedent (□∃x φ(x)):")
    print(f"   - Is ∃x φ(x) true in w1? {antecedent_bf_holds_at_w1}")
    print(f"   - Is ∃x φ(x) true in w2? {antecedent_bf_holds_at_w2}")
    print(f"   => Antecedent is {antecedent_bf_is_true}")

    # Consequent: ∃x □φ(x) at w
    # This means there is an x in D(w) such that for all w' accessible from w,
    # x is in D(w') and φ(x) is true in w'.
    consequent_bf_is_true = False
    for x in D['w']:
        # Check if □φ(x) is true for this x
        x_has_phi_in_w1 = x in D['w1'] and x in I.get('w1', set())
        x_has_phi_in_w2 = x in D['w2'] and x in I.get('w2', set())
        if x_has_phi_in_w1 and x_has_phi_in_w2:
            consequent_bf_is_true = True
            break # Found a witness
    print(f"2. Consequent (∃x □φ(x)):")
    print(f"   - Does 'a' have φ in both w1 and w2? False ('a' is not in D(w2))")
    print(f"   - Does 'b' have φ in both w1 and w2? False ('b' is not in D(w1))")
    print(f"   => Consequent is {consequent_bf_is_true}")

    # Final result for BF
    final_bf_result = not antecedent_bf_is_true or consequent_bf_is_true
    print(f"\nResult for Barcan Formula: {antecedent_bf_is_true} → {consequent_bf_is_true} = {final_bf_result}")
    print("The Barcan formula does not hold in this model.")


    # --- Evaluate Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x) ---
    print("\n\nEvaluating Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)\n")

    # Antecedent of CBF is the same as the consequent of BF
    antecedent_cbf_is_true = consequent_bf_is_true
    print(f"1. Antecedent (∃x □φ(x)): {antecedent_cbf_is_true}")

    # Consequent of CBF is the same as the antecedent of BF
    consequent_cbf_is_true = antecedent_bf_is_true
    print(f"2. Consequent (□∃x φ(x)): {consequent_cbf_is_true}")

    # Final result for CBF
    final_cbf_result = not antecedent_cbf_is_true or consequent_cbf_is_true
    print(f"\nResult for Converse Barcan Formula: {antecedent_cbf_is_true} → {consequent_cbf_is_true} = {final_cbf_result}")
    print("The Converse Barcan formula holds in this model (as False -> True is True).")

if __name__ == '__main__':
    check_formulas()
<<<C>>>