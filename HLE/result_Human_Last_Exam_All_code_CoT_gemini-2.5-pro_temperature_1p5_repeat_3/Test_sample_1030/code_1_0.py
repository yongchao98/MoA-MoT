def p_neg(v):
    """Computes negation in 3-valued logic."""
    mapping = {2: 0, 1: 1, 0: 2} # T:2, G:1, F:0
    return mapping[v]

def p_conj(v1, v2):
    """Computes conjunction (AND) in 3-valued logic."""
    return min(v1, v2)

def p_disj(v1, v2):
    """Computes disjunction (OR) in 3-valued logic."""
    return max(v1, v2)

def p_impl(v1, v2):
    """Computes implication (->) as not v1 or v2."""
    return p_disj(p_neg(v1), v2)

def main():
    """
    Analyzes argument L: A |- (A & B) -> (B & A)
    It demonstrates that the conclusion is a tautology, meaning it's always
    True (2) or Glut (1), and never False (0).
    An argument with a tautological conclusion is always valid.
    """
    print("Analyzing argument L: A |- (A & B) -> (B & A)")
    print("Truth values: T=2 (True), G=1 (Glut), F=0 (False)")
    print("Designated values are {T, G}, i.e., {2, 1}.")
    print("The conclusion is (A & B) -> (B & A). Let's test all possible values for A and B.\n")
    
    vals = {'T': 2, 'G': 1, 'F': 0}
    val_names = {2: 'T', 1: 'G', 0: 'F'}
    is_tautology = True

    for a_name, a_val in vals.items():
        for b_name, b_val in vals.items():
            # Evaluate the antecedent of the conclusion: (A & B)
            p = p_conj(a_val, b_val)
            p_name = val_names[p]

            # The consequent is (B & A), which is identical to (A & B)
            q = p_conj(b_val, a_val)
            q_name = val_names[q]

            # Evaluate the conclusion: P -> Q (which is P -> P here)
            conclusion_val = p_impl(p, q)
            conclusion_name = val_names[conclusion_val]
            
            # Print the step-by-step evaluation
            print(f"Case: A={a_name}, B={b_name}")
            # Output each number in the final equation
            print(f"  (A ∧ B) -> (B ∧ A) becomes ({a_name} ∧ {b_name}) -> ({b_name} ∧ {a_name})")
            print(f"  = ({val_names[a_val]} ∧ {val_names[b_val]}) -> ({val_names[b_val]} ∧ {val_names[a_val]})")
            print(f"  = ({a_val} ∧ {b_val}) -> ({b_val} ∧ {a_val})")
            print(f"  = {p} -> {q}")
            print(f"  = {p_name} -> {q_name}")
            print(f"  = {conclusion_name} (Value: {conclusion_val})")

            # Check if the result is ever False
            if conclusion_val == 0:
                is_tautology = False
                print("  Result is F (False)! This is not a tautology.")
            else:
                print("  Result is designated (T or G).")
            print("-" * 20)

    if is_tautology:
        print("\nConclusion: The formula (A & B) -> (B & A) is a tautology because its value is never F.")
        print("An argument with a tautological conclusion is always valid.")
        print("Therefore, argument L is valid in system KG.")

if __name__ == "__main__":
    main()
