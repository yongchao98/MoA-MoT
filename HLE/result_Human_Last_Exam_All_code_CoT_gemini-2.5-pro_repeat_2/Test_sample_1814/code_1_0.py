import itertools

def powerset(iterable):
    """
    Generates the powerset of an iterable.
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))

def get_saturated_sets(X, T):
    """
    Computes the saturated sets for a topology T on a finite set X.
    A set is saturated if it is the intersection of all its open supersets.
    X: frozenset, the universal set.
    T: set of frozensets, the topology (open sets).
    """
    saturated_sets = set()
    # Iterate over all subsets of X to check which are saturated
    for S in powerset(X):
        S = frozenset(S)
        # Find all open sets in T that are supersets of S
        open_supersets = [U for U in T if S.issubset(U)]
        
        # The intersection of an empty collection of sets is the universe X
        if not open_supersets:
            intersection = X
        else:
            intersection = frozenset.intersection(*open_supersets)
        
        if S == intersection:
            saturated_sets.add(S)
    return saturated_sets

def get_dual_topology(X, T):
    """
    Computes the dual topology d(T) for a topology T on a finite set X.
    On a finite set, compact saturated sets are just saturated sets.
    """
    # 1. The saturated sets form a sub-basis for the closed sets of the dual topology.
    sub_basis_closed = get_saturated_sets(X, T)

    # 2. Generate the basis for closed sets by closing the sub-basis under finite union.
    basis_closed = set(sub_basis_closed)
    while True:
        # Create all possible new unions and check if they are already in the basis
        newly_added = {s1.union(s2) for s1 in basis_closed for s2 in sub_basis_closed} - basis_closed
        if not newly_added:
            break
        basis_closed.update(newly_added)

    # 3. Generate all closed sets by taking all possible intersections of the basis elements.
    closed_sets = set()
    for subset_to_intersect in powerset(basis_closed):
        if not subset_to_intersect:
            # Intersection of an empty collection is the universe
            closed_sets.add(X)
        else:
            closed_sets.add(frozenset.intersection(*subset_to_intersect))

    # 4. The dual topology's open sets are the complements of its closed sets.
    dual_T = {X - closed_set for closed_set in closed_sets}
    
    return dual_T

def print_topology(T, label):
    """Helper function to print a topology in a readable, sorted format."""
    # T is a set of frozensets
    # Convert to a list of sorted lists for consistent, readable output
    printable_sets = sorted([sorted(list(s)) for s in T])
    print(f"{label}: {printable_sets}")

def main():
    """
    Main function to run the dual topology iteration.
    """
    # Define the set X
    X = frozenset({0, 1})

    # Define the initial topology T_0 (Sierpinski space)
    # The open sets are {}, {0}, {0, 1}
    T_0 = {frozenset(), frozenset({0}), frozenset({0, 1})}
    
    # List to store the sequence of distinct topologies found
    topologies = [T_0]
    
    print("Starting with topology T_0:")
    print_topology(T_0, "T_0")
    print("-" * 20)

    current_T = T_0
    i = 1
    while True:
        # Compute the next dual topology in the sequence
        next_T = get_dual_topology(X, current_T)
        
        # Check if this topology has been seen before
        if next_T in topologies:
            print(f"T_{i} is the same as a previously found topology. Stopping.")
            break
        
        # Store and print the new topology
        topologies.append(next_T)
        print_topology(next_T, f"T_{i}")
        
        current_T = next_T
        i += 1
        # Safety break for very long sequences
        if i > 20: 
            print("Stopping after 20 iterations.")
            break
    
    print("-" * 20)
    print(f"The iteration produced a total of {len(topologies)} distinct topologies.")

if __name__ == "__main__":
    main()
