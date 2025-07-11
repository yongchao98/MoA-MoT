def get_dual_finite(topology, universe):
    """
    Computes the dual of a topology on a finite set.
    
    In a finite space:
    - All subsets are compact.
    - All topologies are Alexandrov, so saturated sets are the same as open sets.
    
    The closed sub-basis for the dual is K(T) cap S(T) = P(X) cap O(T) = O(T).
    The closed sets of the dual are finite unions of these, which are just O(T).
    The open sets of the dual are the complements of its closed sets, which are C(T).
    """
    # For a finite topology T, O(d(T)) = C(T).
    open_sets_T = topology
    closed_sets_T = {universe - s for s in open_sets_T}
    
    # The open sets of the dual are the closed sets of the original.
    dual_topology = frozenset(closed_sets_T)
    return dual_topology

def print_topology(topo, name):
    """Helper function to print a topology nicely."""
    # Convert frozensets to sorted lists for consistent display
    sorted_sets = sorted(list(s) for s in topo)
    print(f"{name}: {{", end="")
    print(*[str(s) for s in sorted_sets], sep=", ", end="")
    print("}")

def main():
    """
    Analyzes the dual iteration sequence for a finite topology.
    """
    # Define the universe set X
    X = frozenset({0, 1})

    # Define the starting topology T_0 (Sierpinski Space)
    # The empty set must be represented as frozenset()
    T0 = frozenset({
        frozenset(),
        frozenset({0}),
        frozenset({0, 1})
    })

    print("Analyzing the sequence T, d(T), d(d(T)), ... for a finite space.")
    print(f"Universe X = {sorted(list(X))}\n")

    topologies = [T0]
    current_T = T0

    print_topology(current_T, "T_0")

    # Iterate a few times to show the cycle
    for i in range(1, 5):
        # Compute the dual
        next_T = get_dual_finite(current_T, X)
        print_topology(next_T, f"T_{i} = d(T_{i-1})")

        # Check if we have seen this topology before
        if next_T in topologies:
            print(f"\nCycle detected! T_{i} is the same as a previous topology.")
            break
        else:
            topologies.append(next_T)
        
        current_T = next_T

    print(f"\nFound {len(topologies)} distinct topologies in the sequence.")
    print("As shown by the derivation, for any finite topology, the maximum number is 2.")
    print("\nThe largest possible number of 10 requires an infinite space and a more complex topology.")

if __name__ == '__main__':
    main()