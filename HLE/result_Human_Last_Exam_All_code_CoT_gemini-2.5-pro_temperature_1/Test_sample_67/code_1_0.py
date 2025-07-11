def solve():
    """
    Analyzes the Barcan and Converse Barcan formulas in a Kripke model
    with decreasing domains and prints the findings.
    """
    # 1. Define the Kripke model with decreasing domains
    # Worlds
    W = {'w', 'v1', 'v2'}
    # Accessibility Relation: from 'w', we can access 'v1' and 'v2'
    R = {'w': {'v1', 'v2'}, 'v1': set(), 'v2': set()}
    # Domains: D(v1) and D(v2) are subsets of D(w)
    D = {'w': {'a', 'b'}, 'v1': {'a'}, 'v2': {'b'}}

    # Predicate φ(x) truth values. Assumed False if not specified.
    # phi[world][individual] -> True/False
    phi = {
        'v1': {'a': True},
        'v2': {'b': True}
    }

    # --- Helper functions for evaluating logic ---

    def check_phi(world, individual):
        """Checks if φ(individual) is true at a given world."""
        # An individual must exist in a world to have a property there.
        if individual not in D[world]:
            return False
        return phi.get(world, {}).get(individual, False)

    def check_exists_phi(world):
        """Checks if ∃x φ(x) is true at a world."""
        for individual in D[world]:
            if check_phi(world, individual):
                return True
        return False

    def check_box(world, formula_checker_func):
        """Checks if □(formula) is true at a world."""
        accessible_worlds = R.get(world, set())
        # Vacuously true if no worlds are accessible
        if not accessible_worlds:
            return True
        # Check if the formula holds in ALL accessible worlds
        for accessible_world in accessible_worlds:
            if not formula_checker_func(accessible_world):
                return False
        return True

    # --- Print model details ---
    print("--- Analysis in a Model with Decreasing Domains ---")
    print(f"Worlds: {W}")
    print(f"Domains: {D}")
    print(f"Accessibility from w: R(w) = {R['w']}")
    print(f"Predicate φ truth values: φ(a) is true in v1, φ(b) is true in v2\n")

    # --- Test the Barcan Formula: □∃x φ(x) → ∃x □φ(x) at world 'w' ---
    print("--- Testing the Barcan Formula: □∃x φ(x) → ∃x □φ(x) ---")

    # Evaluate Antecedent: □∃x φ(x) at 'w'
    # This checks if check_exists_phi is true for all worlds accessible from 'w'
    bf_antecedent_holds = check_box('w', check_exists_phi)
    
    print("\n1. Evaluate the antecedent: □∃x φ(x) at 'w'")
    print("   Is it true that in every accessible world (v1, v2), someone has property φ?")
    print(f"   - In 'v1': ∃x φ(x) is {check_exists_phi('v1')} (because of 'a').")
    print(f"   - In 'v2': ∃x φ(x) is {check_exists_phi('v2')} (because of 'b').")
    print(f"   Since it holds for all accessible worlds, the antecedent is: {bf_antecedent_holds}\n")

    # Evaluate Consequent: ∃x □φ(x) at 'w'
    bf_consequent_holds = False
    for individual in D['w']:
        # Check if □φ(individual) holds at 'w'
        if check_box('w', lambda world: check_phi(world, individual)):
            bf_consequent_holds = True
            break
            
    print("2. Evaluate the consequent: ∃x □φ(x) at 'w'")
    print("   Is there a specific individual in D(w) who has property φ in ALL accessible worlds?")
    print(f"   - Check 'a': Does 'a' have φ in v1 AND v2? No, φ(a) is {check_phi('v2', 'a')} in v2.")
    print(f"   - Check 'b': Does 'b' have φ in v1 AND v2? No, φ(b) is {check_phi('v1', 'b')} in v1.")
    print(f"   Since no individual in D(w) satisfies the condition, the consequent is: {bf_consequent_holds}\n")

    # Final result for the Barcan Formula
    bf_holds = not bf_antecedent_holds or bf_consequent_holds
    print("3. Final Result for the Barcan Formula")
    print(f"   The implication is: {bf_antecedent_holds} → {bf_consequent_holds}")
    print(f"   This evaluates to {bf_holds}.")
    print("   Conclusion: The Barcan formula FAILS in this model.\n")

    # --- Test the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x) ---
    print("--- Testing the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x) ---")
    cbf_antecedent_holds = bf_consequent_holds
    cbf_consequent_holds = bf_antecedent_holds
    cbf_holds = not cbf_antecedent_holds or cbf_consequent_holds
    print(f"   The antecedent (∃x □φ(x)) is {cbf_antecedent_holds}.")
    print(f"   The consequent (□∃x φ(x)) is {cbf_consequent_holds}.")
    print(f"   The implication is: {cbf_antecedent_holds} → {cbf_consequent_holds}")
    print(f"   This evaluates to {cbf_holds}.")
    print("   Conclusion: The Converse Barcan Formula HOLDS in this model.\n")

    print("Overall finding: In systems with decreasing domains, the converse Barcan formula holds, but the Barcan formula does not.")

solve()
<<<C>>>