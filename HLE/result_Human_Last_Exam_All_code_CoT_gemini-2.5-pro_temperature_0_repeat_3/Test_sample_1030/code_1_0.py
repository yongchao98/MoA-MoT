import itertools

# Step 1 & 2: Define the 3-valued logic semantics (FDE)
F = frozenset([0])
T = frozenset([1])
G = frozenset([0, 1])
VALS = {'T': T, 'F': F, 'G': G}
VAL_NAMES = {v: k for k, v in VALS.items()}

def is_designated(v):
    """A value is designated if it contains 1 (i.e., is T or G)."""
    return 1 in v

def p_neg(v):
    """Negation flips the presence of 0 and 1."""
    res = set()
    if 0 in v: res.add(1)
    if 1 in v: res.add(0)
    return frozenset(res)

def p_and(v1, v2):
    """Conjunction is true iff both parts are true."""
    res = set()
    if 1 in v1 and 1 in v2: res.add(1)
    if 0 in v1 or 0 in v2: res.add(0)
    return frozenset(res)

def p_or(v1, v2):
    """Disjunction is true if either part is true."""
    res = set()
    if 1 in v1 or 1 in v2: res.add(1)
    if 0 in v1 and 0 in v2: res.add(0)
    return frozenset(res)

def p_implies(v1, v2):
    """Implication A -> B is defined as ¬A ∨ B."""
    return p_or(p_neg(v1), v2)

def check_formula(name, formula_func, num_vars):
    """Checks if a formula is a tautology."""
    print(f"--- Checking Formula {name} ---")
    is_tautology = True
    # Iterate through all possible assignments of T, F, G to variables
    for values in itertools.product(VALS.values(), repeat=num_vars):
        result = formula_func(*values)
        if not is_designated(result):
            is_tautology = False
            var_names = [chr(ord('A') + i) for i in range(num_vars)]
            assignments = {var_names[i]: VAL_NAMES[values[i]] for i in range(num_vars)}
            print(f"Found counter-model: {assignments}")
            print(f"Result: {VAL_NAMES[result]} (Not Designated)\n")
            break
    if is_tautology:
        print(f"Result: Formula {name} is a tautology.\n")
    else:
        print(f"Result: Formula {name} is NOT a tautology.\n")


def check_argument(name, premise_funcs, conclusion_func, num_vars):
    """Checks if an argument is valid."""
    print(f"--- Checking Argument {name} ---")
    is_valid = True
    # Iterate through all possible assignments of T, F, G to variables
    for values in itertools.product(VALS.values(), repeat=num_vars):
        # Check if all premises are designated
        premises_designated = all(is_designated(p_func(*values)) for p_func in premise_funcs)
        
        if premises_designated:
            # If premises are designated, conclusion must also be
            conclusion_val = conclusion_func(*values)
            if not is_designated(conclusion_val):
                is_valid = False
                var_names = [chr(ord('A') + i) for i in range(num_vars)]
                assignments = {var_names[i]: VAL_NAMES[values[i]] for i in range(num_vars)}
                print(f"Found counter-model: {assignments}")
                for i, p_func in enumerate(premise_funcs):
                    print(f"  Premise {i+1} value: {VAL_NAMES[p_func(*values)]} (Designated)")
                print(f"  Conclusion value: {VAL_NAMES[conclusion_val]} (Not Designated)\n")
                break
    if is_valid:
        print(f"Result: Argument {name} is valid.\n")
    else:
        print(f"Result: Argument {name} is NOT valid.\n")

# Step 3 & 4: Define and check the propositional options

# I. ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))
formula_I = lambda vA, vB, vC: p_implies(p_implies(p_or(vA, vB), vC), p_or(p_neg(vA), p_and(p_neg(vB), vC)))
check_formula("I", formula_I, 3)

# G. A → B, B → (¬C ∧ (A ∨ D)) ⊢ A → (¬C ∧ A)
prem_G1 = lambda vA, vB, vC, vD: p_implies(vA, vB)
prem_G2 = lambda vA, vB, vC, vD: p_implies(vB, p_and(p_neg(vC), p_or(vA, vD)))
conc_G = lambda vA, vB, vC, vD: p_implies(vA, p_and(p_neg(vC), vA))
check_argument("G", [prem_G1, prem_G2], conc_G, 4)

# K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
prem_K = lambda vA, vB: p_and(vA, vB)
conc_K = lambda vA, vB: p_implies(p_or(p_neg(vA), p_neg(vB)), p_and(vA, vB))
check_argument("K", [prem_K], conc_K, 2)

# L. A ⊢ (A ∧ B) → (B ∧ A)
prem_L = lambda vA, vB: vA
conc_L = lambda vA, vB: p_implies(p_and(vA, vB), p_and(vB, vA))
check_argument("L", [prem_L], conc_L, 2)

# Step 5: Conclusion based on script output
print("--- Final Conclusion ---")
print("The analysis shows that both arguments K and L are valid in the logic KG.")
print("Argument K, `A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)`, simplifies to `P ⊢ P` (since the conclusion is logically equivalent to the premise). This is the principle of reflexivity, which is the most fundamental property of any logical entailment relation.")
print("Argument L, `A ⊢ (A ∧ B) → (B ∧ A)`, is an instance of an argument where the conclusion is a tautology. The conclusion `(A ∧ B) → (B ∧ A)` is always designated, regardless of the truth values of A and B.")
print("Given that a single choice must be made, and both are technically valid, K represents a more fundamental and universally accepted logical principle (reflexivity of entailment) than L. Therefore, K is the most robustly correct answer.")
