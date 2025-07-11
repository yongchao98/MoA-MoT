def solve_modal_logic_problem():
    """
    Analyzes the Barcan and Converse Barcan formulas in a model
    with decreasing domains to determine their validity.
    """

    # 1. Define the Kripke Model with decreasing domains
    # Worlds: w0 is the starting world, accessible worlds are w1, w2
    # Domains: D(w1) and D(w2) are subsets of D(w0)
    # Predicate phi: Defined to create a counterexample for the Barcan formula
    model = {
        'worlds': {'w0', 'w1', 'w2'},
        'access': {'w0': {'w1', 'w2'}, 'w1': set(), 'w2': set()},
        'domains': {'w0': {'a', 'b'}, 'w1': {'a'}, 'w2': {'b'}},
        'phi': {
            ('a', 'w1'): True, ('a', 'w2'): False,
            ('b', 'w1'): False, ('b', 'w2'): True,
        }
    }

    def get_phi(individual, world):
        return model['phi'].get((individual, world), False)

    print("--- Model with Decreasing Domains ---")
    print(f"Domains: D(w0)={model['domains']['w0']}, D(w1)={model['domains']['w1']}, D(w2)={model['domains']['w2']}")
    print(f"Predicate: phi('a') is True in w1, phi('b') is True in w2\n")

    # 2. Evaluate the Barcan Formula: □∃x φ(x) → ∃x □φ(x) at w0

    print("--- Evaluating Barcan Formula: □∃x φ(x) → ∃x □φ(x) ---")

    # Antecedent: □∃x φ(x)
    # Check if in all worlds accessible from w0, there exists an individual with property phi.
    print("Antecedent (Left Side): □∃x φ(x)")
    # Check w1
    exists_in_w1 = any(get_phi(x, 'w1') for x in model['domains']['w1'])
    print(f"  In accessible world w1, is ∃x φ(x) true? {exists_in_w1}")
    # Check w2
    exists_in_w2 = any(get_phi(x, 'w2') for x in model['domains']['w2'])
    print(f"  In accessible world w2, is ∃x φ(x) true? {exists_in_w2}")
    bf_antecedent = exists_in_w1 and exists_in_w2
    print(f"Result: Antecedent is {bf_antecedent}\n")

    # Consequent: ∃x □φ(x)
    # Check if there exists an individual in w0 that has property phi in all accessible worlds.
    print("Consequent (Right Side): ∃x □φ(x)")
    found_individual = False
    for x in model['domains']['w0']:
        # Check if □φ(x) is true for this x
        phi_in_all_accessible = all(get_phi(x, w_acc) for w_acc in model['access']['w0'])
        print(f"  For individual '{x}' in D(w0), is φ('{x}') true in all accessible worlds (w1, w2)? {phi_in_all_accessible}")
        if phi_in_all_accessible:
            found_individual = True
            break
    bf_consequent = found_individual
    print(f"Result: Consequent is {bf_consequent}\n")

    print("--- Final Evaluation of Barcan Formula ---")
    print(f"The implication is: {bf_antecedent} → {bf_consequent}")
    if bf_antecedent and not bf_consequent:
        print("This is False. The Barcan Formula does NOT hold.\n")
    else:
        print("This is True. The formula holds in this model.\n")


    # 3. Evaluate the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)
    print("--- Evaluating Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x) ---")
    cbf_antecedent = bf_consequent
    cbf_consequent = bf_antecedent
    print(f"Antecedent (Left Side): ∃x □φ(x) is {cbf_antecedent}")
    print(f"Consequent (Right Side): □∃x φ(x) is {cbf_consequent}\n")

    print("--- Final Evaluation of Converse Barcan Formula ---")
    print(f"The implication is: {cbf_antecedent} → {cbf_consequent}")
    if cbf_antecedent and not cbf_consequent:
        print("This is False. The formula does not hold.\n")
    else:
        print("This is True. The Converse Barcan Formula holds in this model (as expected from theory).\n")

    print("--- Overall Conclusion ---")
    print("In systems with decreasing domains:")
    print("The Converse Barcan Formula holds, but the Barcan Formula does not necessarily hold.")

solve_modal_logic_problem()