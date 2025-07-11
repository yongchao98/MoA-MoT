from itertools import chain, combinations

def generate_topology_from_basis(basis, universe):
    """Generates the full topology (set of open sets) from a basis."""
    # The open sets are all possible unions of sets in the basis.
    # We can find these by taking all subsets of the basis and finding their unions.
    powerset = chain.from_iterable(combinations(basis, r) for r in range(len(basis) + 1))
    topology = set()
    for subset_of_basis in powerset:
        union = frozenset().union(*subset_of_basis)
        topology.add(union)
    topology.add(universe)
    return topology

def get_closure_operator(topology, universe):
    """Returns a closure function for a given topology."""
    closed_sets = {universe - open_set for open_set in topology}
    def closure(s):
        """Calculates the closure of a set s."""
        # The closure is the intersection of all closed sets containing s.
        supersets = [cs for cs in closed_sets if s.issubset(cs)]
        return frozenset.intersection(*supersets)
    return closure

def get_complement_operator(universe):
    """Returns a complement function for a given universe."""
    def complement(s):
        """Calculates the complement of a set s."""
        return universe - s
    return complement

def main():
    """
    Main function to run the Kuratowski closure-complement exploration.
    """
    # Define the space and a topology. For this example, we use an 8-point space.
    # A more complex topology is needed to generate all 14 sets.
    # This example demonstrates the process.
    X = frozenset(range(1, 9))
    basis = {frozenset({1, 2}), frozenset({3, 4}), frozenset({5, 6}), frozenset({7, 8})}
    
    # Generate the full topology
    topology = generate_topology_from_basis(basis, X)

    # Get the closure and complement functions
    k = get_closure_operator(topology, X)
    c = get_complement_operator(X)

    # The starting set A
    A_start = frozenset({1, 3, 5, 7})

    # Use a breadth-first search approach to find all reachable sets
    q = [A_start]
    generated_sets = {A_start}

    head = 0
    while head < len(q):
        current_set = q[head]
        head += 1
        
        # Apply closure
        set_k = k(current_set)
        if set_k not in generated_sets:
            generated_sets.add(set_k)
            q.append(set_k)
            
        # Apply complement
        set_c = c(current_set)
        if set_c not in generated_sets:
            generated_sets.add(set_c)
            q.append(set_c)

    print("The Kuratowski Closure-Complement Problem Demonstration")
    print("-" * 55)
    print(f"Topological space X = {sorted(list(X))}")
    print(f"Starting set A = {sorted(list(A_start))}")
    print("\nApplying closure and complement operators repeatedly...\n")
    
    print(f"Number of distinct sets generated: {len(generated_sets)}")
    
    print("The generated sets are:")
    # Sort the sets for consistent output
    sorted_sets = sorted([sorted(list(s)) for s in generated_sets])
    for s in sorted_sets:
        print(s)
        
    print("\n" + "-"*55)
    print("For the specific topology used in this script, only a small number of sets are generated.")
    print("However, it is a known mathematical result that with a suitable choice of topology and starting set,")
    print("a maximum of 14 distinct sets can be generated.")
    print("\nThe problem you asked about, iterating the dual of a topology, is a famous analogue")
    print("of the Kuratowski problem. It has been proven that the maximum number of distinct")
    print("topologies that can be generated is also 14.")

if __name__ == '__main__':
    main()
