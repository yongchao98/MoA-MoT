def print_obstruction_groups():
    """
    This function prints the list of homology and homotopy groups that
    are the fundamental components of the obstruction to homotoping
    the two paths of bundle automorphisms.
    """
    # The variable 'n' represents the dimension of the sphere S^{n-1} which X is a homology-analogue of.
    # The variable 'k' is such that the rank of the vector bundle is 2k.
    n = "n"
    k = "2k"

    # The list of groups that are the building blocks for the obstruction.
    # H_{i}(X) refers to the i-th homology group of X.
    # pi_{j}(SO(2k)) refers to the j-th homotopy group of SO(2k).
    obstruction_groups = [
        f"H_{n}-1(X)",          # From the topology of the base space X
        f"pi_1(SO({k}))",       # From the 'local' automorphisms in the fiber
        f"pi_{n}-1(SO({k}))",   # Related to the classification of the bundle E
        f"pi_{n}(SO({k}))"      # Related to the primary 'global' obstruction
    ]

    print("The obstructions are determined by the following groups:")
    for group in obstruction_groups:
        print(group)

print_obstruction_groups()