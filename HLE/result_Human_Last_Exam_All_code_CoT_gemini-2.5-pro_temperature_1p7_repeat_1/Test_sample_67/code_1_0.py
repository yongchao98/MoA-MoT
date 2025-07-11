def main():
    """
    This script demonstrates that the Barcan formula does not hold in Kripke models
    with decreasing domains by building and evaluating a specific counterexample.
    """
    # 1. Define the Kripke Model with decreasing domains.
    # We have a starting world 'w0' and two accessible worlds, 'w1' and 'w2'.
    worlds = {'w0', 'w1', 'w2'}
    
    # Accessibility Relation: From w0, we can "see" w1 and w2.
    accessible_from = {'w0': ['w1', 'w2'], 'w1': [], 'w2': []}
    
    # Decreasing Domains: The individuals in w1 and w2 are subsets of those in w0.
    domains = {
        'w0': {'a', 'b'},
        'w1': {'a'},       # D(w1) is a subset of D(w0)
        'w2': {'b'}        # D(w2) is a subset of D(w0)
    }
    
    # Valuation for predicate φ: Defines where φ(x) is true.
    # Let φ(a) be true in w1 and φ(b) be true in w2.
    phi_is_true = {('a', 'w1'), ('b', 'w2')}

    print("--- Model Setup for Counterexample ---")
    print(f"Worlds: {worlds}")
    print(f"Domains: {domains}")
    print(f"Accessibility from w0: {accessible_from['w0']}")
    print(f"Truths for φ: φ(a) is true in w1, φ(b) is true in w2\n")

    # --- Helper functions for evaluating formulas ---
    
    def eval_phi(individual, world):
        # φ(x) is true for an individual in a world if the individual
        # exists in that world's domain AND the proposition is valued as true.
        return individual in domains[world] and (individual, world) in phi_is_true

    def eval_exists_phi(world):
        # "There exists an x such that φ(x)" is true in a world if at least one
        # individual in its domain satisfies φ.
        for individual in domains[world]:
            if eval_phi(individual, world):
                return (True, individual) # Return True and the witness
        return (False, None)

    def eval_box(formula_check_func, world_to_check, individual=None):
        # "It is necessary that P" (□P) is true if P is true in all accessible worlds.
        accessible_worlds = accessible_from.get(world_to_check, [])
        if not accessible_worlds: # Vacuously true if no worlds are accessible
            return True
        for next_world in accessible_worlds:
            # Check if the formula fails in any accessible world
            if individual:
                if not formula_check_func(individual, next_world):
                    return (False, next_world) # Fails for the individual in this world
            else:
                if not formula_check_func(next_world)[0]:
                    return (False, next_world) # The existential formula fails in this world
        return (True, None) # Holds in all accessible worlds

    # --- Evaluation of the Barcan Formula at world w0 ---
    print("--- Evaluating Barcan Formula: □∃x φ(x) → ∃x □φ(x) at w0 ---")

    # 1. Evaluate the Antecedent (LHS): □∃x φ(x)
    print("\nStep 1: Evaluate Antecedent (LHS) --> □∃x φ(x)")
    lhs_is_true, failing_world = eval_box(eval_exists_phi, 'w0')
    
    print("  - To be true, ∃x φ(x) must hold in all worlds accessible from w0 (w1, w2).")
    exists_w1, witness_w1 = eval_exists_phi('w1')
    print(f"  - Checking w1: Does ∃x φ(x) hold? {'Yes' if exists_w1 else 'No'}. (Witness: {witness_w1})")
    exists_w2, witness_w2 = eval_exists_phi('w2')
    print(f"  - Checking w2: Does ∃x φ(x) hold? {'Yes' if exists_w2 else 'No'}. (Witness: {witness_w2})")
    
    if lhs_is_true:
        print("  --> LHS is TRUE because ∃x φ(x) holds in all accessible worlds.")
    else:
        print(f"  --> LHS is FALSE because ∃x φ(x) fails in world {failing_world}.")

    # 2. Evaluate the Consequent (RHS): ∃x □φ(x)
    print("\nStep 2: Evaluate Consequent (RHS) --> ∃x □φ(x)")
    print("  - To be true, we need to find an individual in D(w0) = {a, b} for whom □φ(x) is true.")
    
    rhs_is_true = False
    witness_for_rhs = None

    for individual in domains['w0']:
        print(f"  - Checking individual '{individual}': Is □φ({individual}) true at w0?")
        print(f"    - This means φ({individual}) must be true in all worlds accessible from w0 (w1, w2).")
        is_necessary, failing_world = eval_box(eval_phi, 'w0', individual=individual)
        if is_necessary:
            print(f"    --> YES, □φ({individual}) is TRUE.")
            rhs_is_true = True
            witness_for_rhs = individual
            break
        else:
            print(f"    --> NO, □φ({individual}) is FALSE because φ({individual}) fails in world {failing_world}.")

    if rhs_is_true:
         print(f"  --> RHS is TRUE. (Witness: {witness_for_rhs})")
    else:
         print("  --> RHS is FALSE because no individual from D(w0) satisfies φ in all accessible worlds.")
    
    # 3. Final Result for the implication
    print("\n--- Final Result ---")
    final_value = (not lhs_is_true) or rhs_is_true
    print(f"The Barcan Formula is the implication: {lhs_is_true} → {rhs_is_true}")
    print(f"This implication is {'TRUE' if final_value else 'FALSE'} in our model.")
    print("\nBecause we found a counterexample, the Barcan formula does not generally hold for decreasing domains.")

if __name__ == '__main__':
    main()
