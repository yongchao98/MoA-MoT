# A simple Kripke model representation for first-order modal logic
class KripkeModel:
    def __init__(self, worlds, access_relation, domains, assignments):
        self.worlds = worlds
        self.access_relation = access_relation
        self.domains = domains
        self.assignments = assignments # assignments for predicate 'phi'

    def get_accessible_worlds(self, world):
        return self.access_relation.get(world, [])

    # Check if a property 'phi' is true for an individual in a world
    def phi(self, world, individual):
        return self.assignments.get(world, {}).get(individual, False)

# Functions to evaluate modal formulas
def evaluate_exists_phi(model, world):
    """Evaluates ∃x φ(x) in a given world."""
    for individual in model.domains.get(world, []):
        if model.phi(world, individual):
            return True, [individual] # Return true and the witness
    return False, []

def evaluate_box_phi(model, world, individual):
    """Evaluates □φ(x) for a given individual in a world."""
    if individual not in model.domains.get(world, []):
        # The individual must exist in the starting world to even ask this
        return False

    for accessible_world in model.get_accessible_worlds(world):
        # The individual must exist in the accessible world and phi must be true
        if individual not in model.domains.get(accessible_world, []) or \
           not model.phi(accessible_world, individual):
            return False
    return True

def evaluate_exists_box_phi(model, world):
    """Evaluates ∃x □φ(x) in a given world."""
    for individual in model.domains.get(world, []):
        if evaluate_box_phi(model, world, individual):
            return True, [individual] # Return true and the witness
    return False, []

def evaluate_box_exists_phi(model, world):
    """Evaluates □∃x φ(x) in a given world."""
    witnesses = {}
    for accessible_world in model.get_accessible_worlds(world):
        result, witness = evaluate_exists_phi(model, accessible_world)
        if not result:
            return False, {} # A world was found where nothing has phi
        witnesses[accessible_world] = witness
    return True, witnesses


# --- Main Program ---
if __name__ == "__main__":
    # Define the Kripke model (counter-example for Barcan Formula)
    # Worlds: w0, w1, w2
    # Accessibility: w0 -> w1, w0 -> w2
    # Domains (decreasing): D(w0) = {a, b}, D(w1) = {a}, D(w2) = {b}
    # Predicate 'phi' assignments: φ(a) is true in w1, φ(b) is true in w2
    model = KripkeModel(
        worlds={'w0', 'w1', 'w2'},
        access_relation={'w0': ['w1', 'w2']},
        domains={'w0': ['a', 'b'], 'w1': ['a'], 'w2': ['b']},
        assignments={'w1': {'a': True}, 'w2': {'b': True}}
    )

    start_world = 'w0'
    print("--- Analysis of Formulas at World w0 ---")
    print(f"Model has decreasing domains: D(w1) ⊆ D(w0) and D(w2) ⊆ D(w0)\n")

    # 1. Evaluate Barcan Formula: □∃x φ(x) → ∃x □φ(x)
    print("1. Barcan Formula: □∃x φ(x) → ∃x □φ(x)")

    # Antecedent: □∃x φ(x)
    ant_bf_val, ant_bf_witnesses = evaluate_box_exists_phi(model, start_world)
    print(f"   - Antecedent (□∃x φ(x)) is: {ant_bf_val}")
    print(f"     (Because in w1, ∃x φ(x) is true for x={ant_bf_witnesses['w1'][0]}, and in w2, ∃x φ(x) is true for x={ant_bf_witnesses['w2'][0]})")


    # Consequent: ∃x □φ(x)
    con_bf_val, con_bf_witnesses = evaluate_exists_box_phi(model, start_world)
    print(f"   - Consequent (∃x □φ(x)) is: {con_bf_val}")
    print(f"     (Because no single individual from {model.domains['w0']} exists and has property φ in all accessible worlds)")

    # Full implication
    holds_bf = not ant_bf_val or con_bf_val
    print(f"   => Result: The Barcan Formula is {holds_bf} for this model.\n")

    # 2. Evaluate Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)
    print("2. Converse Barcan Formula: ∃x □φ(x) → □∃x φ(x)")

    # Antecedent: ∃x □φ(x) (This is the same as the consequent from before)
    ant_cbf_val = con_bf_val
    print(f"   - Antecedent (∃x □φ(x)) is: {ant_cbf_val}")


    # Consequent: □∃x φ(x) (This is the same as the antecedent from before)
    con_cbf_val = ant_bf_val
    print(f"   - Consequent (□∃x φ(x)) is: {con_cbf_val}")


    # Full implication
    holds_cbf = not ant_cbf_val or con_cbf_val
    print(f"   => Result: The Converse Barcan Formula is {holds_cbf} for this model.")

    print("\n--- Conclusion ---")
    print("The simulation provides a counter-example showing the Barcan Formula can fail in a system with decreasing domains.")
    print("The Converse Barcan Formula holds (as the antecedent is false, the implication is vacuously true in this specific case, and holds generally by logical proof).")
    print("Therefore, the Converse Barcan formula holds, but the Barcan formula does not hold in all possible worlds.")
