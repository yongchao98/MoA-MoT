def main():
    """
    Determines the truth value of the given statement by constructing a counterexample.
    """
    
    # Let's define the three-valued logic implication.
    # We use a standard definition where p -> q is 1 if p <= q, and 0 if p=1, q<1.
    # And 0.5 otherwise (e.g. p=0.5, q=0)
    def implies(p, q):
        if p <= q:
            return 1.0
        elif p == 1.0 and q < 1.0:
            return 0.0
        elif p == 0.5 and q == 0.0:
            return 0.5
        # This covers all cases for {0, 0.5, 1}
        return 1.0

    print("Step 1: Define the problem.")
    print("We need to find the truth value of S: Box(P) in world w1, where")
    print("P = forall x, y, z: (T(x, y, z) -> Box(T(x, y, z)))")
    print("-" * 30)

    print("Step 2: Find a counterexample for the inner formula P.")
    print("Axiom ATV only fixes the truth value of T(x, y, z) across worlds when 'z' is a world like w1, w2, etc.")
    print("It does NOT constrain T(x, y, z) if 'z' is a non-world individual like 'a'.")
    print("Let's choose (x, y, z) = (1, 'b', 'a') and define a model consistent with the axioms.")
    print("-" * 30)
    
    # The set of worlds in the equivalence class.
    worlds = ['w1', 'w2', 'w3']
    
    # Our counterexample assignment for the predicate T(1, 'b', 'a').
    # This is consistent with all axioms because 'a' is a concrete entity, not a world.
    model_T_1_b_a = {
        'w1': 1.0,
        'w2': 0.0,
        'w3': 1.0
    }

    print("Step 3: Evaluate the implication T(1, b, a) -> Box(T(1, b, a)) in world w1.")
    
    # Antecedent: T(1, b, a) evaluated in w1
    antecedent_val = model_T_1_b_a['w1']
    print(f"The truth value of the antecedent T(1, 'b', 'a') in w1 is: {antecedent_val}")
    
    # Consequent: Box(T(1, b, a)) evaluated in w1
    # This is the minimum value of T(1, 'b', 'a') across all accessible worlds.
    consequent_val = min(model_T_1_b_a.values())
    print(f"The truth value of the consequent Box(T(1, 'b', 'a')) in w1 is min({model_T_1_b_a['w1']}, {model_T_1_b_a['w2']}, {model_T_1_b_a['w3']}) = {consequent_val}")
    
    # Evaluate the implication
    implication_val = implies(antecedent_val, consequent_val)
    print(f"The truth value of the implication ({antecedent_val} -> {consequent_val}) is: {implication_val}")
    print("-" * 30)
    
    print("Step 4: Determine the truth value of P in w1.")
    # Since we found a case where the implication is false (0), the universally quantified
    # statement P is false in w1.
    P_val_in_w1 = implication_val # In a universal, one false instance makes the whole thing false.
    print(f"Since a counterexample exists, the value of P in w1 is {P_val_in_w1}.")
    print("-" * 30)

    print("Step 5: Determine the truth value of Box(P) in w1.")
    # The value of Box(P) is the minimum value of P across all accessible worlds.
    # Since P is 0 in w1, the minimum can be at most 0.
    final_truth_value = P_val_in_w1 # min(P_in_w1, P_in_w2, ...) = min(0, ...) = 0
    print(f"The value of Box(P) in w1 is min(Value(P) in w1, Value(P) in w2, ...) = min({P_val_in_w1}, ...)")
    print(f"Therefore, the final truth value is {final_truth_value}.")

if __name__ == "__main__":
    main()