def get_obstruction_groups(n, k):
    """
    Identifies the homotopy-theoretic obstruction groups.

    Args:
        n: The dimension of the homology sphere X is (n-1).
        k: The rank of the vector bundle E is 2k.
    """
    
    # The obstruction to homotopy between the two paths lies in the group
    # pi_1(Aut(E)), where Aut(E) is the space of bundle automorphisms of E.
    # This group is isomorphic to the product of two groups.

    # Group 1: This component arises from the global rotation of the fibers.
    # It corresponds to the fundamental group of the structure group SO(2k).
    group1 = f"pi_1(SO({2*k}))"

    # Group 2: This component depends on the topology of the base space Sigma(X).
    # It is given by [Sigma^2(X), SO(2k)]_*.
    # Since X is a homology (n-1)-sphere, Sigma^2(X) is homotopy equivalent
    # to an (n+1)-sphere, S^(n+1).
    # Thus, the group is the (n+1)-th homotopy group of SO(2k).
    group2 = f"pi_{n+1}(SO({2*k}))"

    print("The homotopy-theoretic obstructions are classified by elements in two groups.")
    print("These groups are:")
    print(f"1. {group1}")
    print(f"2. {group2}")

# Example usage:
# Let's say X is a homology 3-sphere (n-1=3 => n=4)
# and the bundle E is rank 8 (2k=8 => k=4).
n_example = 4
k_example = 4
get_obstruction_groups(n_example, k_example)
