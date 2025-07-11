import itertools

def get_saturated_sets(topology, X):
    """
    Computes the collection of saturated sets for a given topology.
    A set is saturated if it is an intersection of open sets.
    On a finite space, compact sets are all subsets, so we only need to
    care about saturated sets for the dual operation.
    """
    # The collection of open sets is `topology`.
    # We will find all possible intersections of these open sets.
    saturated_sets = set()
    
    # An intersection of an empty collection of sets is the whole space X.
    saturated_sets.add(X)

    topo_list = list(topology)
    # Consider all non-empty sub-collections of the topology
    for r in range(1, len(topo_list) + 1):
        for sub_collection in itertools.combinations(topo_list, r):
            # Calculate the intersection of the open sets in this sub-collection
            # Start with the whole space and intersect down.
            intersection = set(X)
            for open_set in sub_collection:
                intersection.intersection_update(open_set)
            saturated_sets.add(frozenset(intersection))
            
    return saturated_sets

def compute_dual_topology(topology, X):
    """
    Computes the dual of a given topology.
    The dual topology has a closed sub-basis consisting of all
    compact saturated sets of the original topology.
    On a finite space, this simplifies:
    - All subsets are compact.
    - So, compact saturated sets are just saturated sets.
    The dual topology has an OPEN sub-basis consisting of the
    complements of the saturated sets.
    """
    # 1. Find all saturated sets.
    saturated_sets = get_saturated_sets(topology, X)

    # 2. Form the open sub-basis for the dual topology.
    # This is the collection of complements of the saturated sets.
    open_sub_basis = {frozenset(X - s) for s in saturated_sets}

    # 3. Generate the open basis from the sub-basis.
    # The basis is the set of all finite intersections of sub-basis elements.
    open_basis = set()
    open_basis.add(X) # The empty intersection is X
    sub_basis_list = list(open_sub_basis)
    for r in range(1, len(sub_basis_list) + 1):
        for combo in itertools.combinations(sub_basis_list, r):
            intersection = set(X)
            for s in combo:
                intersection.intersection_update(s)
            open_basis.add(frozenset(intersection))

    # 4. Generate the full dual topology from the basis.
    # The topology is the set of all arbitrary unions of basis elements.
    dual_topology = set()
    basis_list = list(open_basis)
    for r in range(len(basis_list) + 1):
        for combo in itertools.combinations(basis_list, r):
            union = set()
            for s in combo:
                union.update(s)
            dual_topology.add(frozenset(union))
            
    return dual_topology

def main():
    """
    Main function to run the demonstration.
    """
    print("Exploring the iteration of the dual topology on a finite space.")
    print("This illustrates the concept, but the maximum number (7) is achieved on an infinite space.\n")

    # Define a finite set X
    X = frozenset({0, 1, 2})

    # Define an initial topology T0 on X. A topology is a set of frozensets.
    # T0 must contain the empty set and X, and be closed under union and finite intersection.
    T0 = {
        frozenset(),
        frozenset({0}),
        frozenset({1}),
        frozenset({0, 1}),
        frozenset({0, 1, 2})
    }
    
    # Store the sequence of generated topologies
    topologies = [T0]
    current_topology = T0

    # Iterate the dual operation until a cycle is found
    for i in range(1, 10): # Limit iterations to prevent infinite loops in case of error
        next_topology = compute_dual_topology(current_topology, X)
        
        print(f"T{i-1} has {len(current_topology)} open sets.")
        # To make output readable, we can print the sets themselves for small examples
        # print(sorted([sorted(list(s)) for s in current_topology]))

        if next_topology in topologies:
            start_index = topologies.index(next_topology)
            print(f"\nT{i} is the same as T{start_index}.")
            print("\nA cycle has been found.")
            # The final equation here is showing which topologies are equal
            print(f"The sequence of distinct topologies is [T0, T1, ..., T{len(topologies)-1}].")
            print(f"The number of distinct topologies found for this example is {len(topologies)}.")
            print(f"The sequence becomes periodic: ..., T{start_index}, T{start_index+1}, ...")
            return
        
        topologies.append(next_topology)
        current_topology = next_topology
        
if __name__ == '__main__':
    main()