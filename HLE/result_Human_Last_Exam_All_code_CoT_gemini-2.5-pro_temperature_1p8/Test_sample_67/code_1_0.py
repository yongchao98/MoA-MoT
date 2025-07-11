import itertools

class KripkeModel:
    """A simple Kripke model for first-order modal logic."""
    def __init__(self, worlds, relation, domains, interpretation):
        self.W = worlds  # Set of worlds
        self.R = relation  # Accessibility relation: dict mapping world -> set of accessible worlds
        self.D = domains  # Domain function: dict mapping world -> set of objects
        self.I = interpretation # Interpretation function: dict mapping (property, world) -> set of objects

    def check_phi(self, obj, world):
        """Checks if phi(obj) is true at world."""
        # Standard "serious actualist" semantics: for a property to be true of an object,
        # the object must exist in that world's domain.
        if obj not in self.D[world]:
            return False
        return obj in self.I.get((phi, world), set())

# Let's define the property 'phi' for clarity in the interpretation dict key
def phi(x):
    return f"phi({x})"

# ---- 1. Setup of the Counterexample Model for the Barcan Formula ----
# This model has decreasing domains: D(w1) and D(w2) are subsets of D(w0)
print("Setting up a Kripke model with decreasing domains...\n")
model = KripkeModel(
    worlds={'w0', 'w1', 'w2'},
    relation={'w0': {'w1', 'w2'}, 'w1': set(), 'w2': set()},
    domains={'w0': {'a', 'b'}, 'w1': {'a'}, 'w2': {'b'}},
    interpretation={
        (phi, 'w1'): {'a'}, # In world w1, 'a' has property phi
        (phi, 'w2'): {'b'}  # In world w2, 'b' has property phi
    }
)

# ---- 2. Evaluate the Barcan Formula: □∃x φ(x) → ∃x □φ(x) at w0 ----
print("--- Evaluating the Barcan Formula (BF): □∃x φ(x) → ∃x □φ(x) at w0 ---")

# Evaluate Antecedent: □∃x φ(x) at w0
# "In every world accessible from w0, does there exist an object with property phi?"
antecedent_holds = True
print("Evaluating Antecedent: □∃x φ(x)")
accessible_worlds = model.R['w0']
if not accessible_worlds:
    antecedent_holds = True # Vacuously true if no worlds are accessible
    print("  No worlds accessible from w0, so antecedent is vacuously true.")
else:
    for world in accessible_worlds:
        # Check if ∃x φ(x) is true in 'world'
        exists_phi_in_world = any(model.check_phi(obj, world) for obj in model.D[world])
        print(f"  Checking accessible world '{world}': Is ∃x φ(x) true?")
        print(f"    Domain D({world}) = {model.D[world]}")
        if exists_phi_in_world:
            witness = next(obj for obj in model.D[world] if model.check_phi(obj, world))
            print(f"    Yes, object '{witness}' is a witness.")
        else:
            print(f"    No, no object in D({world}) satisfies φ.")
            antecedent_holds = False
            break

print(f"\nResult for Antecedent (□∃x φ(x)): {antecedent_holds}\n")


# Evaluate Consequent: ∃x □φ(x) at w0
# "Does there exist an object in w0 that has property phi in ALL accessible worlds?"
consequent_holds = False
print("Evaluating Consequent: ∃x □φ(x)")
print(f"  Domain D(w0) = {model.D['w0']}")
print("  We need to find an object in D(w0) that has property φ in ALL accessible worlds.")
for obj in model.D['w0']:
    print(f"  Testing candidate object '{obj}' from D(w0)...")
    # Check if □φ(obj) is true at w0 for this obj
    is_necessary_for_obj = True
    for world in accessible_worlds:
        if model.check_phi(obj, world):
            print(f"    - In world '{world}', '{obj}' has property φ. (OK)")
        else:
            reason = f"object '{obj}' does not exist in D({world})." if obj not in model.D[world] else f"I(φ, {world}) does not include '{obj}'."
            print(f"    - In world '{world}', '{obj}' does NOT have property φ because {reason} (FAIL)")
            is_necessary_for_obj = False
            break
    if is_necessary_for_obj:
        print(f"  SUCCESS: Object '{obj}' works for all accessible worlds.")
        consequent_holds = True
        break
    else:
        print(f"  FAILURE: Object '{obj}' does not work for all accessible worlds.")

print(f"\nResult for Consequent (∃x □φ(x)): {consequent_holds}\n")

# Final check of the implication for BF
bf_is_true = not antecedent_holds or consequent_holds
print(f"Barcan Formula evaluates to: {antecedent_holds} → {consequent_holds}, which is {bf_is_true}.")
if not bf_is_true:
    print("This counterexample shows the Barcan Formula does NOT hold in systems with decreasing domains.")

# ---- 3. Evaluate the Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x) ----
print("\n--- Evaluating the Converse Barcan Formula (CBF): ∃x □φ(x) → □∃x φ(x) at w0 ---")
# Using the same model. The antecedent '∃x □φ(x)' was already shown to be False.
# In logic, if the antecedent is false, the implication (False -> X) is always true.
print("From our evaluation above, we know the antecedent '∃x □φ(x)' is FALSE.")
print("The consequent '□∃x φ(x)' is TRUE.")
cbf_is_true = not consequent_holds or antecedent_holds # Using the values from the BF evaluation
cbf_is_true_rearranged = not consequent_holds or antecedent_holds
cbf_is_true = not consequent_holds or antecedent_holds # Antecedent of CBF is consequent of BF
antecedent_cbf = consequent_holds
consequent_cbf = antecedent_holds

cbf_is_true = not antecedent_cbf or consequent_cbf
print(f"Converse Barcan Formula evaluates to: {antecedent_cbf} → {consequent_cbf}, which is {cbf_is_true}.")
print("The Converse Barcan Formula holds in this model. In fact, it can be proven to hold in all models with decreasing domains.")
print("\nConclusion:")
print("In systems with decreasing domains:")
print("- The Barcan formula does not hold.")
print("- The converse Barcan formula holds.")
<<<C>>>