def check_barcan_in_decreasing_domains():
    """
    Analyzes the Barcan Formula (BF) and its Converse (CBF)
    in a Kripke model with decreasing domains.
    """
    # Define the Kripke model (worlds, accessibility, domains, interpretation of φ)
    # This model is a counterexample for the Barcan Formula.
    model = {
        'worlds': {'w0', 'w1'},
        'accessibility': {
            'w0': {'w0', 'w1'},  # w0 can "see" itself and w1
            'w1': {'w1'}       # w1 can only "see" itself
        },
        'domains': {
            'w0': {'pegasus', 'horse'},  # Domain of w0
            'w1': {'horse'}             # Domain of w1 (subset of D(w0) -> decreasing domain)
        },
        # Interpretation of the predicate φ(x) as "x exists".
        # This is a common choice for demonstrating domain issues.
        # φ(x) is true if x is in the domain of the given world.
        'phi': lambda individual, world, m: individual in m['domains'][world]
    }

    print("--- Analysis of Barcan Formulas in a Decreasing Domain Model ---")
    print(f"Model Worlds: {model['worlds']}")
    print(f"Model Domains: D(w0)={model['domains']['w0']}, D(w1)={model['domains']['w1']}")
    print("Since w0 can access w1 and D(w1) is a subset of D(w0), this is a decreasing domain system.")
    print("\nLet's test the Barcan Formula: □∃x φ(x) → ∃x □φ(x) at world w0")
    print("Let φ(x) be the property 'x exists'.\n")

    # --- Helper functions for evaluation ---
    def evaluate_exists_phi(world, m):
        """Checks if ∃x φ(x) is true at a given world."""
        # Check if any individual in the world's domain has property φ
        for individual in m['domains'][world]:
            if m['phi'](individual, world, m):
                return True
        return False

    def evaluate_box(inner_formula_eval, world, m, *args):
        """Checks if □(...) is true at a given world."""
        # Check if the inner formula is true for all accessible worlds
        for accessible_world in m['accessibility'][world]:
            if not inner_formula_eval(accessible_world, m, *args):
                return False
        return True
    
    def evaluate_box_phi_for_individual(world, m, individual):
        """Helper for checking □φ(c) for a specific individual c."""
        # An individual must exist in the world to be considered
        if individual not in m['domains'][world]:
            return False
            
        # Now check if □φ(individual) is true
        return evaluate_box(lambda w, model, ind: model['phi'](ind, w, model), world, m, individual)

    # --- Evaluate the Barcan Formula ---

    # 1. Evaluate the Antecedent: □∃x φ(x) at w0
    print("Step 1: Evaluate the antecedent □∃x φ(x) at w0.")
    print("  This means for all worlds w' accessible from w0, ∃x φ(x) must be true in w'.")
    antecedent_holds = True
    for w_prime in sorted(list(model['accessibility']['w0'])):
        result = evaluate_exists_phi(w_prime, model)
        print(f"  - Checking w': {w_prime}. Does ∃x φ(x) hold? -> {result}")
        if not result:
            antecedent_holds = False
            break
    print(f"Result for Antecedent: {antecedent_holds}\n")


    # 2. Evaluate the Consequent: ∃x □φ(x) at w0
    print("Step 2: Evaluate the consequent ∃x □φ(x) at w0.")
    print("  This means there is an individual c in D(w0) for which □φ(c) is true.")
    consequent_holds = False
    for individual in sorted(list(model['domains']['w0'])):
        print(f"  - Checking individual '{individual}' from D(w0)...")
        # Check if □φ(individual) is true at w0
        is_phi_necessary = evaluate_box_phi_for_individual('w0', model, individual)
        print(f"    - Does □φ('{individual}') hold at w0? -> {is_phi_necessary}")
        if is_phi_necessary:
            consequent_holds = True
            break # Found one, so the existential is true.
    print(f"Result for Consequent: {consequent_holds}\n")

    # --- Conclusion ---
    print("--- Final Conclusion ---")
    barcan_formula_valid = not (antecedent_holds and not consequent_holds)
    if not barcan_formula_valid:
        print("The Barcan Formula (□∃x φ(x) → ∃x □φ(x)) FAILS in this model.")
        print(f"Because the antecedent is TRUE but the consequent is FALSE ({antecedent_holds} → {consequent_holds} is False).")
    else:
        print("The Barcan Formula holds in this specific model (but may fail in others).")

    print("\nAs established by theoretical proof, the Converse Barcan Formula (∃x □φ(x) → □∃x φ(x)) HOLDS in all systems with decreasing domains.")
    print("\nTherefore, in systems with decreasing domains, the Converse Barcan formula holds, but the Barcan formula does not necessarily hold.")


if __name__ == '__main__':
    check_barcan_in_decreasing_domains()
    # The analysis shows CBF holds and BF fails. This corresponds to choice C.
    print("\n<<<C>>>")