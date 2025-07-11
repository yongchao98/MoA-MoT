# Define truth values
T = 2  # True
G = 1  # Glut (Both True and False)
F = 0  # False
vals = {'T': T, 'G': G, 'F': F}
val_names = {v: k for k, v in vals.items()}
DESIGNATED = {T, G}

# Define logical connectives for the 3-valued logic KG
def negation(v):
    """ ¬v """
    return 2 - v

def conjunction(a, b):
    """ a ∧ b """
    return min(a, b)

def disjunction(a, b):
    """ a ∨ b """
    return max(a, b)

def implication(a, b):
    """ a → b, defined as ¬a ∨ b """
    return disjunction(negation(a), b)

def check_argument_L():
    """
    Checks the validity of argument L: A ⊢ (A ∧ B) → (B ∧ A)
    An argument is valid if it's impossible for the premise to be designated
    and the conclusion to be non-designated.
    """
    print("Checking validity of argument L: A ⊢ (A ∧ B) → (B ∧ A)")
    print("-" * 50)
    print("v(A)\tv(B)\tPremise: A\tConclusion: (A ∧ B) → (B ∧ A)\tValid?")
    
    is_valid = True
    
    # Iterate over all possible truth values for atomic propositions A and B
    for v_A in vals.values():
        for v_B in vals.values():
            
            # Evaluate the premise
            premise_val = v_A
            
            # Evaluate the conclusion
            # Let X = A ∧ B and Y = B ∧ A. The conclusion is X → Y
            # Since conjunction is commutative, v(X) will equal v(Y)
            val_A_and_B = conjunction(v_A, v_B)
            val_B_and_A = conjunction(v_B, v_A)
            conclusion_val = implication(val_A_and_B, val_B_and_A)

            # Check for invalidity condition
            premise_is_designated = premise_val in DESIGNATED
            conclusion_is_designated = conclusion_val in DESIGNATED
            
            status = "OK"
            if premise_is_designated and not conclusion_is_designated:
                is_valid = False
                status = "INVALID"
            
            print(f"{val_names[v_A]}\t{val_names[v_B]}\t"
                  f"{val_names[premise_val]} ({'des' if premise_is_designated else 'undes'})\t\t"
                  f"{val_names[conclusion_val]} ({'des' if conclusion_is_designated else 'undes'})\t\t\t{status}")

    print("-" * 50)
    if is_valid:
        print("Result: The argument is VALID.")
        print("The conclusion (A ∧ B) → (B ∧ A) is a tautology, as its value is always T or G (designated).")
        print("An argument with a tautologous conclusion is always valid.")
    else:
        print("Result: The argument is INVALID.")

check_argument_L()