def check_barcan_formulas():
    """
    This function models a Kripke structure with decreasing domains to test
    the Barcan formula and its converse.
    """
    # 1. Define the Kripke Model M = <W, R, D, I>
    # W = Worlds
    worlds = {'w1', 'w2', 'w3'}
    # R = Accessibility Relation
    # From w1, we can "see" w2 and w3.
    # We only care about accessibility from w1 for this example.
    relation = {'w1': ['w2', 'w3']}
    # D = Domain function (Decreasing Domains)
    # D(w2) is a subset of D(w1), D(w3) is a subset of D(w1).
    domains = {
        'w1': {'a', 'b'},
        'w2': {'a'},
        'w3': {'b'},
    }
    # I = Interpretation function for a predicate P, which stands for φ
    # We define where P(x) is true.
    predicate_p = {
        'w1': set(),
        'w2': {'a'}, # P(a) is true in w2
        'w3': {'b'}, # P(b) is true in w3
    }

    # Helper function to check if a formula 'phi(x)' holds for an individual 'x' in a 'world'
    def holds(world, x):
        return x in predicate_p[world]

    # Helper function to evaluate Necessarily (□)
    def necessarily(world, formula_check_func):
        if world not in relation:
            return True # Vacuously true if no worlds are accessible
        for accessible_world in relation[world]:
            if not formula_check_func(accessible_world):
                return False
        return True

    # Helper function to evaluate Exists (∃)
    def exists(world, formula_check_func):
        for individual in domains[world]:
            if formula_check_func(world, individual):
                return True
        return False

    print("--- Analysis of the Barcan Formula in a Decreasing Domain Model ---")
    print(f"Model: W={worlds}, D(w1)={domains['w1']}, D(w2)={domains['w2']}, D(w3)={domains['w3']}")
    print("Accessibility: w1 -> w2, w1 -> w3")
    print("Property φ (as P): P(a) is true in w2. P(b) is true in w3.")
    print("\nBarcan Formula: □∃x φ(x) → ∃x □φ(x)")

    # 2. Evaluate the antecedent: □∃x φ(x) at w1
    # Define the inner formula: ∃x φ(x)
    def exists_phi(world):
        return exists(world, holds)

    antecedent_holds = necessarily('w1', exists_phi)

    print(f"\n1. Checking the Antecedent (Left Side): □∃x φ(x) at w1")
    print("   Is it true that in all worlds accessible from w1, there exists an x such that φ(x)?")
    print(f"   - In w2: ∃x φ(x)? Yes, 'a' exists and has property φ. Value: {exists_phi('w2')}")
    print(f"   - In w3: ∃x φ(x)? Yes, 'b' exists and has property φ. Value: {exists_phi('w3')}")
    print(f"   Since it holds for all accessible worlds, the antecedent is: {antecedent_holds}")


    # 3. Evaluate the consequent: ∃x □φ(x) at w1
    def exists_nec_phi(world):
        # For each individual 'c' in the current world's domain...
        for c in domains[world]:
            # ...check if it necessarily has property phi.
            # Define a check for □φ(c)
            def nec_phi_c(accessible_world):
                # Under actualism, if c doesn't exist in the world, the property is false for it.
                if c not in domains[accessible_world]:
                    return False
                return holds(accessible_world, c)

            if necessarily(world, nec_phi_c):
                # If we find even one such individual, the existential is true.
                return True, c
        return False, None

    consequent_holds, witness = exists_nec_phi('w1')

    print(f"\n2. Checking the Consequent (Right Side): ∃x □φ(x) at w1")
    print("   Is there an individual in w1 (a or b) that necessarily has property φ?")
    print("   - Check 'a': Is □φ(a) true? Does 'a' have property φ in w2 AND w3?")
    print(f"     - φ(a) in w2? {holds('w2', 'a')}")
    print(f"     - φ(a) in w3? {'a' not in domains['w3']} -> False")
    print("     So, □φ(a) is False.")
    print("   - Check 'b': Is □φ(b) true? Does 'b' have property φ in w2 AND w3?")
    print(f"     - φ(b) in w2? {'b' not in domains['w2']} -> False")
    print(f"     - φ(b) in w3? {holds('w3', 'b')}")
    print("     So, □φ(b) is False.")
    print(f"   Since no individual in w1 has the property necessarily, the consequent is: {consequent_holds}")


    # 4. Final conclusion
    print("\n--- Conclusion for Barcan Formula ---")
    print(f"The formula is an implication: {antecedent_holds} → {consequent_holds}")
    if antecedent_holds and not consequent_holds:
        print("This evaluates to TRUE → FALSE, which is FALSE.")
        print("Therefore, the Barcan formula does not hold in this model.")
    else:
        print("This counterexample did not falsify the formula.")


    # 5. Brief check for Converse Barcan Formula
    # CBF: ∃x □φ(x) → □∃x φ(x)
    # The antecedent, ∃x □φ(x), is already false in our model.
    # An implication with a false antecedent (FALSE → ...) is always true.
    # The Converse Barcan Formula holds (as proven in the thinking steps).

    print("\n--- Conclusion for the Problem ---")
    print("The Converse Barcan Formula holds in decreasing domain models.")
    print("The Barcan Formula fails in decreasing domain models.")
    print("The correct answer choice reflects this.")

check_barcan_formulas()