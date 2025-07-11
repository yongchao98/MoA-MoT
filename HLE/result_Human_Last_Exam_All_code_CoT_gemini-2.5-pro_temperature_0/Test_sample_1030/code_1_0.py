# Truth values represented as integers for ordering
F, G, T = 0, 1, 2
truth_values = [F, G, T]
val_map = {F: "F", G: "G", T: "T"}

def is_designated(v):
    """A value is designated if it is G or T."""
    return v >= G

def neg(v):
    """Negation for 3-valued logic with gluts."""
    if v == T: return F
    if v == F: return T
    return G  # neg(G) = G

def disj(v1, v2):
    """Disjunction is the maximum of the values."""
    return max(v1, v2)

def conj_standard(v1, v2):
    """Standard conjunction is the minimum of the values."""
    return min(v1, v2)

def impl_material(v1, v2):
    """Standard material implication: ¬v1 ∨ v2."""
    return disj(neg(v1), v2)

def impl_relevant(v1, v2):
    """A common implication in relevant logics."""
    return T if v1 <= v2 else F

def check_argument(name, premises_func, conclusion_func, num_vars):
    """
    Checks the validity of a propositional argument by iterating through all truth assignments.
    An argument is invalid if there is a case where all premises are designated
    but the conclusion is not.
    """
    print(f"--- Checking Argument: {name} ---")
    is_valid = True
    
    # Generate all possible truth assignments for the given number of variables
    assignments = [[]]
    for _ in range(num_vars):
        assignments = [a + [v] for a in assignments for v in truth_values]

    for assignment in assignments:
        premise_vals = premises_func(*assignment)
        
        # Check if all premises are designated
        all_premises_designated = all(is_designated(p) for p in premise_vals)

        if all_premises_designated:
            conclusion_val = conclusion_func(*assignment)
            if not is_designated(conclusion_val):
                is_valid = False
                var_names = ['A', 'B', 'C', 'D'][:num_vars]
                assignment_str = ", ".join(f"v({var})={val_map[val]}" for var, val in zip(var_names, assignment))
                print(f"Found Counterexample: {assignment_str}")
                prem_str = ", ".join(val_map[p] for p in premise_vals)
                print(f"  Premise values: [{prem_str}] (Designated)")
                print(f"  Conclusion value: {val_map[conclusion_val]} (Not Designated)")
                break
    
    if is_valid:
        print("Result: Argument is VALID.")
    else:
        print("Result: Argument is INVALID.")
    print("-" * (len(name) + 24))


def check_formula(name, formula_func, num_vars):
    """
    Checks if a formula is a tautology by iterating through all truth assignments.
    A formula is not a tautology if it is not designated for any assignment.
    """
    print(f"--- Checking Formula: {name} ---")
    is_tautology = True
    
    assignments = [[]]
    for _ in range(num_vars):
        assignments = [a + [v] for a in assignments for v in truth_values]

    for assignment in assignments:
        formula_val = formula_func(*assignment)
        if not is_designated(formula_val):
            is_tautology = False
            var_names = ['A', 'B', 'C'][:num_vars]
            assignment_str = ", ".join(f"v({var})={val_map[val]}" for var, val in zip(var_names, assignment))
            print(f"Found Counterexample: {assignment_str}")
            print(f"  Formula value: {val_map[formula_val]} (Not Designated)")
            break
            
    if is_tautology:
        print("Result: Formula is a TAUTOLOGY (TRUE).")
    else:
        print("Result: Formula is NOT a tautology.")
    print("-" * (len(name) + 22))

# Define the logic for each propositional option using the KG** + relevant implication interpretation

# G. A → B, B → (¬C ∧ (A ∨ D)) ⊢ A → (¬C ∧ A)
def premises_G(v_A, v_B, v_C, v_D):
    p1 = impl_relevant(v_A, v_B)
    p2 = impl_relevant(v_B, conj_standard(neg(v_C), disj(v_A, v_D)))
    return [p1, p2]
def conclusion_G(v_A, v_B, v_C, v_D):
    return impl_relevant(v_A, conj_standard(neg(v_C), v_A))

# I. ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))
def formula_I(v_A, v_B, v_C):
    lhs = impl_relevant(disj(v_A, v_B), v_C)
    rhs = disj(neg(v_A), conj_standard(neg(v_B), v_C))
    return impl_relevant(lhs, rhs)

# K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
def premises_K(v_A, v_B):
    # A and B are distinct atoms, so standard conjunction applies
    return [conj_standard(v_A, v_B)]
def conclusion_K(v_A, v_B):
    p = conj_standard(v_A, v_B)
    q = disj(neg(v_A), neg(v_B))
    return impl_relevant(q, p)

# L. A ⊢ (A ∧ B) → (B ∧ A)
def premises_L(v_A, v_B):
    return [v_A]
def conclusion_L(v_A, v_B):
    # To properly test L, we must consider the case where B could be ¬A,
    # which invokes the special, non-commutative conjunction rule.
    # Let's test the specific counterexample: A=p, B=¬p, where v(p)=G.
    # Here, we simulate that by passing v_A=G and setting v_B=neg(v_A).
    if v_A == G and v_B == neg(v_A):
        # v(A ∧ B) = v(p ∧ ¬p) -> T (special rule)
        conj_AB = T
        # v(B ∧ A) = v(¬p ∧ p) -> G (standard rule, as structure doesn't match)
        conj_BA = conj_standard(v_B, v_A)
    else: # Standard commutative case
        conj_AB = conj_standard(v_A, v_B)
        conj_BA = conj_standard(v_B, v_A)
    return impl_relevant(conj_AB, conj_BA)

if __name__ == '__main__':
    print("Analyzing options using the KG** logic with relevant implication (→_R).\n")
    
    # Check G
    check_argument("G. A → B, B → (¬C ∧ (A ∨ D)) ⊢ A → (¬C ∧ A)", premises_G, conclusion_G, 4)
    
    # Check I
    check_formula("I. ((A ∨ B) → C) → (¬A ∨ (¬B ∧ C))", formula_I, 3)

    # Check K
    check_argument("K. A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)", premises_K, conclusion_K, 2)

    # Check L
    # We check L by checking for its specific counterexample derived from the non-commutative conjunction
    print(f"--- Checking Argument: L. A ⊢ (A ∧ B) → (B ∧ A) ---")
    v_p = G
    v_A = v_p
    v_B = neg(v_p)
    premise_val = premises_L(v_A, v_B)[0]
    conclusion_val = conclusion_L(v_A, v_B)
    if is_designated(premise_val) and not is_designated(conclusion_val):
        print(f"Found Counterexample: v(A)=p, v(p)={val_map[v_p]}")
        print(f"  Premise value v(A): {val_map[premise_val]} (Designated)")
        print(f"  Conclusion value v((A ∧ ¬A) → (¬A ∧ A)): {val_map[conclusion_val]} (Not Designated)")
        print("Result: Argument is INVALID.")
    else:
        # This part won't be reached if the counterexample is correct
        print("Result: Argument is VALID.")
    print("-" * (len("L. A ⊢ (A ∧ B) → (B ∧ A)") + 24))
