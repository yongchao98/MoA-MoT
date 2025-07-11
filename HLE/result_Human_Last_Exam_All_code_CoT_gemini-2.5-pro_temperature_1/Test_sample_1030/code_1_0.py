def solve():
    """
    This script verifies the validity of argument K in the 3-valued logic KG,
    interpreted as the Logic of Paradox (LP) with an entailment relation based on
    the preservation of truth-degree (v(Premise) <= v(Conclusion)).
    """
    
    # Mapping numeric values to their symbolic representation
    val_map = {0: 'F', 1: 'G', 2: 'T'}

    print("Testing validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print("Using LP semantics and v(Premise) <= v(Conclusion) entailment.\n")
    print("Truth Values (numeric): F=0, G=1, T=2")
    print("Truth Values (symbolic): F, G, T with order F < G < T")
    print("Negation (¬): ¬T=F, ¬G=G, ¬F=T")
    print("Implication (A → B): ¬A ∨ B\n")

    is_valid = True

    # Iterate through all 3x3=9 possible truth value combinations for A and B
    for v_a in range(3):
        for v_b in range(3):
            
            # Define the LP connectives based on numeric values
            def neg(v):
                if v == 2: return 0  # ¬T = F
                if v == 1: return 1  # ¬G = G
                if v == 0: return 2  # ¬F = T

            def conj(v1, v2):
                return min(v1, v2)

            def disj(v1, v2):
                return max(v1, v2)

            def impl(v1, v2):
                return disj(neg(v1), v2)

            # --- Evaluation for the current case ---
            
            # 1. Calculate the value of the Premise: A ∧ B
            v_premise = conj(v_a, v_b)

            # 2. Calculate the value of the Conclusion: (¬A ∨ ¬B) → (A ∧ B)
            # Part 1 of conclusion: (¬A ∨ ¬B)
            v_neg_a = neg(v_a)
            v_neg_b = neg(v_b)
            v_disj_negs = disj(v_neg_a, v_neg_b)

            # Part 2 of conclusion: (A ∧ B) - same as premise
            v_conj_ab = v_premise
            
            # Final conclusion value
            v_conclusion = impl(v_disj_negs, v_conj_ab)

            # --- Output the evaluation steps ---
            
            print(f"Case: v(A)={val_map[v_a]} ({v_a}), v(B)={val_map[v_b]} ({v_b})")
            
            # Premise equation
            print(f"  Premise: v(A ∧ B) = v({val_map[v_a]} ∧ {val_map[v_b]}) = {v_premise}")

            # Conclusion equation
            print(f"  Conclusion: v((¬A ∨ ¬B) → (A ∧ B))")
            print(f"    v(¬{val_map[v_a]} ∨ ¬{val_map[v_b]}) = v({val_map[v_neg_a]} ∨ {val_map[v_neg_b]}) = {v_disj_negs}")
            print(f"    v({val_map[v_disj_negs]} → {val_map[v_conj_ab]}) = v(¬{val_map[v_disj_negs]} ∨ {val_map[v_conj_ab]}) = v({val_map[neg(v_disj_negs)]} ∨ {val_map[v_conj_ab]}) = {v_conclusion}")
            
            # 3. Check the validity condition for this case
            check_result = v_premise <= v_conclusion
            print(f"  Validity Check: v(Premise) <= v(Conclusion)?")
            print(f"    {v_premise} <= {v_conclusion}  ->  {check_result}")
            print("-" * 30)

            if not check_result:
                is_valid = False

    if is_valid:
        print("\nResult: The argument is VALID.")
        print("The condition v(Premise) <= v(Conclusion) holds true for all possible value assignments.")
    else:
        print("\nResult: The argument is INVALID.")
        print("A counterexample was found where v(Premise) > v(Conclusion).")

solve()