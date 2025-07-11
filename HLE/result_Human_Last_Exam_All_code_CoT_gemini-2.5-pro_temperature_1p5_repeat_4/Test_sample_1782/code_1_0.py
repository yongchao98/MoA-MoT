import itertools

def find_branches(tower, path, level, branches):
    """Recursively finds all branches in a tower of partitions."""
    if level == len(tower):
        branches.append(path)
        return

    # The current element of the path we are extending
    parent_node = path[-1]
    
    # Iterate through elements in the current level's partition
    for child_node in tower[level]:
        # A child node extends the path if it's a subset of the parent
        if child_node.issubset(parent_node):
            find_branches(tower, path + [child_node], level + 1, branches)

def main():
    """
    Demonstrates that a tower of refining partitions has a common refinement
    constructed from the infima of its branches.
    """
    base_set = frozenset(range(4))
    print(f"Base set X = {set(base_set)}\n")

    # Create a tower of partitions (maximal antichains) L_0, L_1, L_2
    # where L_{i+1} refines L_i
    l0 = {frozenset(range(4))}
    l1 = {frozenset({0, 1}), frozenset({2, 3})}
    l2 = {frozenset({0}), frozenset({1}), frozenset({2, 3})}
    l3 = {frozenset({0}), frozenset({1}), frozenset({2}), frozenset({3})}
    
    tower = [l0, l1, l2, l3]

    print("Tower of partitions:")
    for i, level in enumerate(tower):
        print(f"L_{i}: {[set(s) for s in level]}")

    # Find all branches through the tower
    branches = []
    # Start the search from the root of the tree (the single element in L0)
    find_branches(tower, [list(tower[0])[0]], 1, branches)
    
    print("\nFound branches (sequences x_i in L_i where x_{i+1} subset x_i):")
    for i, branch in enumerate(branches):
        print(f"  Branch {i+1}: {[set(s) for s in branch]}")

    # Compute the infimum (intersection) of each branch
    infima = {frozenset.intersection(*branch) for branch in branches}
    
    # Filter out empty sets, though none are expected here
    common_refinement = {inf for inf in infima if inf}

    print("\nConstructed common refinement (set of infima of branches):")
    # This is the "final equation" part
    print_friendly_refinement = sorted([set(s) for s in common_refinement], key=lambda x: min(x))
    print(f"L_common = {print_friendly_refinement}")

    # Verify that it is a maximal antichain
    print("\nVerifying that L_common is a maximal antichain:")
    
    # Check for pairwise disjointness
    is_disjoint = True
    for s1, s2 in itertools.combinations(common_refinement, 2):
        if not s1.isdisjoint(s2):
            is_disjoint = False
            break
    
    if is_disjoint:
        print("- The elements are pairwise disjoint.")
    else:
        print("- The elements are NOT pairwise disjoint.")

    # Check if the union is the whole set
    union_of_refinement = frozenset.union(*common_refinement)
    
    if union_of_refinement == base_set:
        print(f"- The union of elements is {set(union_of_refinement)}, which equals the base set X.")
        # Another "final equation"
        print(f"  Equation: Union({', '.join(map(str, print_friendly_refinement))}) = {set(base_set)}")

    else:
        print(f"- The union of elements is {set(union_of_refinement)}, which does NOT equal the base set X.")

    print("\nThis illustrates how a common refinement is guaranteed by the structure of the tree.")

if __name__ == '__main__':
    main()