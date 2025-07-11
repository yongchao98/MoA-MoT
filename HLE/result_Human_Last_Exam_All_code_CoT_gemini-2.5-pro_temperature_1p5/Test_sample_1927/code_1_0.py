import collections

class KripkeModel:
    """Represents a Kripke model for intuitionistic logic."""

    def __init__(self, worlds, relations, valuation):
        """
        Initializes the Kripke model.
        - worlds: A set of world names (e.g., {'w0', 'w1', ...}).
        - relations: A dict where keys are world names and values are sets of
                     worlds directly accessible (e.g., {'w0': {'w1'}}).
        - valuation: A dict where keys are prop. variables and values are sets
                     of worlds where the variable is true.
        """
        self.worlds = worlds
        self.relations = collections.defaultdict(set, relations)
        self.valuation = collections.defaultdict(set, valuation)
        self._successors_cache = {}

    def get_successors(self, world):
        """
        Returns all worlds accessible from a given world (reflexive and
        transitive closure).
        """
        if world in self._successors_cache:
            return self._successors_cache[world]

        visited = set()
        queue = collections.deque([world])
        visited.add(world)
        while queue:
            current = queue.popleft()
            for successor in self.relations.get(current, []):
                if successor not in visited:
                    visited.add(successor)
                    queue.append(successor)
        self._successors_cache[world] = visited
        return visited

    def forces(self, world, formula):
        """
        Checks if a world forces a given formula.
        The formula is represented as a tuple.
        """
        # Base case for 'False'
        if formula == 'False':
            return False
        # Case for propositional variables
        if isinstance(formula, str):
            return world in self.valuation[formula]

        op = formula[0]
        # Conjunction (AND)
        if op == 'and':
            return self.forces(world, formula[1]) and self.forces(world, formula[2])
        # Disjunction (OR)
        if op == 'or':
            return self.forces(world, formula[1]) or self.forces(world, formula[2])
        # Implication (->)
        if op == '->':
            sub_formula1, sub_formula2 = formula[1], formula[2]
            return all(not self.forces(w, sub_formula1) or self.forces(w, sub_formula2)
                       for w in self.get_successors(world))
        raise ValueError(f"Unknown formula structure: {formula}")

def construct_countermodel():
    """Constructs the 8-node Kripke countermodel."""
    # Define the structure of the 8-node model
    worlds = {f'w{i}' for i in range(8)}
    relations = {
        'w0': {'w1'},
        'w1': {'w2', 'w3'},
        'w2': {'w4', 'w5'},
        'w3': {'w6', 'w7'}
    }
    
    # Define the valuation based on the derivation
    # V(A0) = {w4, w6}
    # V(B0) = {}
    # V(A1) = {w2, w4, w5} (to satisfy monotonicity)
    # V(B1) = {w4, w5, w6, w7}
    # V(B2) = {w2, w3, w4, w5, w6, w7} (to satisfy w'|=P1 => w'|=B2 for w'>=w1)
    valuation = {
        'A0': {'w4', 'w6'},
        'B0': set(),
        'A1': {'w2', 'w4', 'w5'},
        'B1': {'w4', 'w5', 'w6', 'w7'},
        'B2': {'w2', 'w3', 'w4', 'w5', 'w6', 'w7'},
    }
    
    return KripkeModel(worlds, relations, valuation)

def main():
    """
    Main function to construct the model, verify it, and print the result.
    """
    # Define the formula structure
    # P0 = (A0 -> B0) v (~A0 -> B0)
    # ~A0 is (A0 -> False)
    Not_A0 = ('->', 'A0', 'False')
    P0 = ('or', ('->', 'A0', 'B0'), ('->', Not_A0, 'B0'))
    
    # P1 = (A1 -> B1) v (~A1 -> B1)
    Not_A1 = ('->', 'A1', 'False')
    P1 = ('or', ('->', 'A1', 'B1'), ('->', Not_A1, 'B1'))

    # Antecedent = [ P0 -> B1 ] ^ [ P1 -> B2 ]
    Antecedent = ('and', ('->', P0, 'B1'), ('->', P1, 'B2'))
    
    # Full formula = Antecedent -> B2
    Formula = ('->', Antecedent, 'B2')

    model = construct_countermodel()
    
    # Verify that w0 does NOT force the formula, making it a countermodel
    if not model.forces('w0', Formula):
        num_nodes = len(model.worlds)
        # To show why it's a countermodel, we show that w1 forces the antecedent but not B2
        w1_forces_antecedent = model.forces('w1', Antecedent)
        w1_forces_B2 = model.forces('w1', 'B2')
        
        # print(f"Verification: The model has {num_nodes} nodes.")
        # print(f"Formula is forced at w0: {model.forces('w0', Formula)}")
        # print(f"w1 forces Antecedent: {w1_forces_antecedent}")
        # print(f"w1 forces B2: {w1_forces_B2}")
        
        # Print the final answer as required
        print(num_nodes)
    else:
        print("Error: The constructed model is not a countermodel.")

if __name__ == "__main__":
    main()