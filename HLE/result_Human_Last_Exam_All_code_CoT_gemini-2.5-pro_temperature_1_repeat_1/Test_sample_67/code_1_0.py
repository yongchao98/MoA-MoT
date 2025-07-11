# A simple Kripke model to test the Barcan formulas.

# --- Model Definition ---

# Worlds
worlds = {'w', 'w_prime'}

# Accessibility Relation (w can 'see' w_prime)
# We assume reflexivity for □, so a world can always see itself.
accessibility = {
    'w': {'w', 'w_prime'},
    'w_prime': {'w_prime'}
}

# Domains for each world (Decreasing: D(w_prime) is a subset of D(w))
domains = {
    'w': {'Alice', 'Bob'},
    'w_prime': {'Alice'}
}

# Predicate extension for phi (e.g., phi(x) means "x is happy")
# (individual, world) -> True/False
phi_extension = {
    ('Alice', 'w'): False,
    ('Bob', 'w'): True,
    ('Alice', 'w_prime'): True,
    # ('Bob', 'w_prime') is implicitly False as Bob is not in D(w_prime)
}

def is_phi_true(individual, world):
    """Checks if phi(individual) is true at a given world."""
    if individual not in domains[world]:
        return False
    return phi_extension.get((individual, world), False)

# --- Formula Evaluation Functions ---

def eval_exists_phi(world):
    """Evaluates ∃x φ(x) at a world."""
    for individual in domains[world]:
        if is_phi_true(individual, world):
            # Found a witness
            return True
    return False

def eval_box(world, formula_func):
    """Evaluates □(formula) at a world."""
    for accessible_world in accessibility[world]:
        if not formula_func(accessible_world):
            return False
    return True

def eval_box_phi_for_individual(individual, world):
    """Evaluates □φ(c) for a specific individual c, from the perspective of a world."""
    for accessible_world in accessibility[world]:
        if not is_phi_true(individual, accessible_world):
            return False
    return True

def eval_exists_box_phi(world):
    """Evaluates ∃x □φ(x) at a world."""
    for individual in domains[world]:
        if eval_box_phi_for_individual(individual, world):
            # Found a witness
            return True
    return False

# --- Main Logic: Test Formulas at World 'w' ---

print("--- Analysis in a Model with Decreasing Domains ---\n")

print("Model setup:")
print(f"Worlds: {worlds}")
print(f"Domains: {domains}")
print(f"Accessibility from 'w': {accessibility['w']}")
print("Predicate 'phi':")
print(f"  - phi(Alice) at 'w': {is_phi_true('Alice', 'w')}")
print(f"  - phi(Bob) at 'w': {is_phi_true('Bob', 'w')}")
print(f"  - phi(Alice) at 'w_prime': {is_phi_true('Alice', 'w_prime')}")
print("-" * 20)

# 1. Test the Barcan Formula: □∃x φ(x) → ∃x □φ(x)

print("\n1. Testing the Barcan Formula: □∃x φ(x) → ∃x □φ(x)")

# Antecedent: □∃x φ(x)
print("Evaluating antecedent: □∃x φ(x) at 'w'")
# Check ∃x φ(x) at all worlds accessible from 'w'
val_exists_phi_at_w = eval_exists_phi('w')
print(f"  - Is ∃x φ(x) true at 'w'? {val_exists_phi_at_w} (Witness: Bob)")
val_exists_phi_at_w_prime = eval_exists_phi('w_prime')
print(f"  - Is ∃x φ(x) true at 'w_prime'? {val_exists_phi_at_w_prime} (Witness: Alice)")
bf_antecedent_holds = eval_box('w', eval_exists_phi)
print(f"Result: Antecedent □∃x φ(x) is {bf_antecedent_holds}\n")

# Consequent: ∃x □φ(x)
print("Evaluating consequent: ∃x □φ(x) at 'w'")
# Check if any individual in D(w) has □φ property
val_box_phi_for_alice = eval_box_phi_for_individual('Alice', 'w')
print(f"  - Is □φ(Alice) true? {val_box_phi_for_alice} (because phi(Alice) is false at 'w')")
val_box_phi_for_bob = eval_box_phi_for_individual('Bob', 'w')
print(f"  - Is □φ(Bob) true? {val_box_phi_for_bob} (because Bob doesn't exist at 'w_prime')")
bf_consequent_holds = eval_exists_box_phi('w')
print(f"Result: Consequent ∃x □φ(x) is {bf_consequent_holds}\n")

# Final result for Barcan Formula
bf_holds = (not bf_antecedent_holds) or bf_consequent_holds
print(f"Barcan Formula holds? {bf_antecedent_holds} → {bf_consequent_holds} is {bf_holds}.")
print("Conclusion: The Barcan Formula is falsified in this model.\n")
print("-" * 20)


# 2. Test the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)
# The reasoning for this is general and doesn't need a specific model check, but we can do it.
# In our model, the antecedent is false, so the implication is vacuously true.
# The general proof in the text above shows it always holds.

print("\n2. Testing the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)")
cbf_antecedent_holds = eval_exists_box_phi('w')
print(f"Evaluating antecedent: ∃x □φ(x) at 'w'. Result is {cbf_antecedent_holds}")
cbf_consequent_holds = eval_box('w', eval_exists_phi)
print(f"Evaluating consequent: □∃x φ(x) at 'w'. Result is {cbf_consequent_holds}")
cbf_holds = (not cbf_antecedent_holds) or cbf_consequent_holds
print(f"Converse Barcan Formula holds? {cbf_antecedent_holds} → {cbf_consequent_holds} is {cbf_holds}.")
print("Conclusion: The Converse Barcan Formula holds (as it always does).")
