def get_obstruction_groups():
    """
    Identifies and prints the list of groups that form the homotopy-theoretic
    obstructions for the described problem.

    The obstructions to the homotopy between the two paths lie in the group
    pi_1 of the space of bundle automorphisms. This group is filtered, and its
    filtration quotients (the actual obstruction groups) are constructed from
    a set of fundamental groups derived from the topology of the space X and
    the special orthogonal group SO(2k).
    """

    # The fundamental groups from which the obstructions are constructed.
    # The notation uses LaTeX format for clarity.
    group_1 = "H_{n-1}(X)"
    group_2 = "\\pi_1(SO(2k))"
    group_3 = "\\pi_{n+1}(SO(2k))"

    print("The homotopy-theoretic obstructions are elements of groups which are constructed from the following list of fundamental groups:")
    print(f"1. The (n-1)-th homology group of X: ${group_1}$")
    print(f"2. The fundamental group of the special orthogonal group SO(2k): ${group_2}$")
    print(f"3. The (n+1)-th homotopy group of the special orthogonal group SO(2k): ${group_3}$")

if __name__ == '__main__':
    get_obstruction_groups()