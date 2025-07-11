def find_nonabelian_filled_groups():
    """
    This script presents the solution to the group theory problem:
    "What are the nonabelian filled groups of order 2q^m for odd prime q and natural number m?"
    The solution is derived from established theorems in group theory, specifically the classification
    of 'filled' or 'C-groups'.
    """

    print("The nonabelian filled groups of order 2 * q^m, where q is an odd prime and m is a natural number, are identified as follows:")
    print("-" * 80)
    print("Based on a classification theorem for filled groups (also known as C-groups), a nonabelian group has this property if and only if it possesses a unique minimal normal subgroup N, and the quotient G/N is a cyclic group of prime order.")
    print("\nAnalysis of groups with order 2 * q^m shows that for a group to meet these criteria, the following conditions must be met:")
    print("1. The parameter 'm' must be equal to 1.")
    print("2. The group must be the semidirect product of a cyclic group of order q and a cyclic group of order 2, with a nontrivial action.")
    print("\nThis uniquely identifies the family of groups up to isomorphism.")
    print("-" * 80)
    print("\nCONCLUSION:")
    print("The only nonabelian filled groups of order 2 * q^m are the DIHEDRAL GROUPS of order 2q. These exist only when m = 1.")
    print("\nFor any odd prime q:")
    print(f"The group is D_2q, and its order is given by the equation: 2 * q^1 = {2}*q.")

    print("\nExamples:")
    q_example_1 = 3
    order_1 = 2 * q_example_1**1
    print(f"For q = {q_example_1}, m = 1: The group is the Dihedral group D_{2*q_example_1} (D_6), with order 2 * {q_example_1}^1 = {order_1}. This group is also isomorphic to the symmetric group S_3.")

    q_example_2 = 5
    order_2 = 2 * q_example_2**1
    print(f"For q = {q_example_2}, m = 1: The group is the Dihedral group D_{2*q_example_2} (D_10), with order 2 * {q_example_2}^1 = {order_2}.")
    
    q_example_3 = 13
    order_3 = 2 * q_example_3**1
    print(f"For q = {q_example_3}, m = 1: The group is the Dihedral group D_{2*q_example_3} (D_26), with order 2 * {q_example_3}^1 = {order_3}.")
    
    print("\nIf m > 1, no nonabelian filled groups of order 2 * q^m exist.")

if __name__ == '__main__':
    find_nonabelian_filled_groups()