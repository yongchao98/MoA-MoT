from itertools import chain, combinations

def powerset(iterable):
    """
    Computes the powerset of an iterable.
    Example: powerset([1,2]) --> () (1,) (2,) (1,2)
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def get_saturated_sets(topology, universe):
    """
    Calculates the saturated sets for a given topology.
    A set is saturated if it is an intersection of open sets.
    """
    open_sets = list(topology)
    saturated_sets = {universe} # The intersection of an empty collection of sets is the universe
    
    # Calculate all possible non-empty intersections of open sets
    for subset_of_opens in powerset(open_sets):
        if not subset_of_opens:
            continue
        
        # Calculate the intersection of the current subset of open sets
        intersection = universe.copy()
        for open_set in subset_of_opens:
            intersection.intersection_update(open_set)
        saturated_sets.add(frozenset(intersection))
        
    return saturated_sets

def dual_operator(topology, universe):
    """
    Applies the dual operator to a topology on a finite set.
    """
    # On a finite space, all subsets are compact.
    # So, the compact saturated sets are just the saturated sets.
    compact_saturated_sets = get_saturated_sets(topology, universe)
    
    # This forms the closed sub-basis for the new topology.
    closed_sub_basis = compact_saturated_sets
    
    # Generate the closed basis: all possible unions of sets from the sub-basis.
    # On a finite space, "finite union" means we can take any number of unions.
    closed_basis = set()
    for subset_of_sub_basis_sets in powerset(closed_sub_basis):
        union = frozenset()
        for s in subset_of_sub_basis_sets:
            union = union.union(s)
        closed_basis.add(union)

    # Generate the closed sets: all possible intersections of sets from the basis.
    closed_sets = set()
    for subset_of_basis_sets in powerset(closed_basis):
        if not subset_of_basis_sets:
            intersection = frozenset(universe) # Intersection of empty collection is the universe
        else:
            # Start with a copy of the first set in the tuple
            iter_subset = iter(subset_of_basis_sets)
            intersection = next(iter_subset).copy()
            # Intersect with the rest
            for s in iter_subset:
                intersection.intersection_update(s)
        closed_sets.add(frozenset(intersection))

    # The new topology is the set of complements of the closed sets.
    new_topology = {frozenset(universe - c) for c in closed_sets}
    
    return new_topology

def format_topology(t, universe_size):
    """Helper function to print topologies in a readable format."""
    return sorted([sorted(list(s)) for s in t])

def main():
    """
    Main function to run the simulation.
    """
    # Define the universe X = {0, 1}
    X = frozenset({0, 1})
    
    # Start with the Sierpinski topology on X
    # T_0 = { {}, {0}, {0,1} }
    T_0 = {frozenset(), frozenset({0}), frozenset({0, 1})}
    
    print(f"Starting with topology T_0 = {format_topology(T_0, len(X))}\n")
    
    distinct_topologies = []
    current_topology = T_0
    
    for i in range(10): # Limit iterations to prevent infinite loops in case of error
        if current_topology in distinct_topologies:
            break
        distinct_topologies.append(current_topology)
        
        print(f"T_{i} = {format_topology(current_topology, len(X))}")
        
        current_topology = dual_operator(current_topology, X)

    print("\n--- Results ---")
    print("The sequence of generated topologies becomes periodic.")
    print(f"The number of distinct topologies found in this sequence is: {len(distinct_topologies)}")

if __name__ == "__main__":
    main()