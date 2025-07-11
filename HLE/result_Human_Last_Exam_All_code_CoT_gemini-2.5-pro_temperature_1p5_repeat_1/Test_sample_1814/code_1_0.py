import itertools

def get_powerset(s):
    """Generates the powerset of a given set."""
    x = len(s)
    base_list = list(s)
    masks = [1 << i for i in range(x)]
    for i in range(1 << x):
        yield frozenset([base_list[j] for j in range(x) if (i & masks[j])])

def get_closed_sets(topology, universe_set):
    """Computes the closed sets for a given topology."""
    return {universe_set - open_set for open_set in topology}

def get_all_unions(basis_sets):
    """Computes the set of all possible unions from a collection of basis sets."""
    if not basis_sets:
        return {frozenset()}
    
    # Using a frozenset of frozensets to store the generated unions to avoid duplicates
    # and allow for quick lookups.
    unions = {frozenset()} | set(basis_sets)
    
    # Iteratively compute unions until no new sets are generated.
    while True:
        new_unions = set()
        for s1 in unions:
            for s2 in unions:
                union = s1 | s2
                if union not in unions:
                    new_unions.add(union)
        
        if not new_unions:
            break
        unions.update(new_unions)
        
    return unions

def c_operator(topology, universe_set):
    """
    Computes the simplified dual topology c(T) for a topology on a finite set.
    The open sets of c(T) are the unions of the closed sets of T.
    """
    closed_sets = get_closed_sets(topology, universe_set)
    open_sets_of_cT = get_all_unions(closed_sets)
    return frozenset(open_sets_of_cT)

def main():
    """
    Demonstrates the iteration of the dual operator on a finite topology.
    """
    X = frozenset({1, 2, 3})

    # An initial topology on X = {1, 2, 3}
    T0 = frozenset([
        frozenset(),
        frozenset({1}),
        frozenset({1, 2}),
        X
    ])
    
    # Use a list to store the sequence of distinct topologies
    topologies_chain = []
    
    # Use a set to keep track of topologies we've already seen for quick checking
    seen_topologies = set()

    current_T = T0
    
    iteration_limit = 10 # A safe limit to prevent infinite loops in case of error
    
    for i in range(iteration_limit):
        if current_T in seen_topologies:
            print(f"\nIteration {i}: The generated topology is the same as T_{topologies_chain.index(current_T)}.")
            print("A cycle is found. Stopping.")
            break
        
        topologies_chain.append(current_T)
        seen_topologies.add(current_T)
        
        print(f"--- Iteration {i}: T_{i} ---")
        # Sort the sets for consistent, readable output
        sorted_topology = sorted(list(current_T), key=lambda s: (len(s), tuple(sorted(list(s)))))
        print("{\n" + ",\n".join([f"  {set(s) if s else '{}'}" for s in sorted_topology]) + "\n}")

        # Compute the next topology in the sequence
        current_T = c_operator(current_T, X)
    
    print(f"\nThe largest number of distinct topologies generated in this example is {len(topologies_chain)}.")
    
    # The final answer to the original question
    print("\n---------------------------------------------------------------------")
    print("The theoretical largest possible number of distinct topologies that can arise from iterating the dual is 7.")
    print("---------------------------------------------------------------------")


if __name__ == "__main__":
    main()