def find_nonabelian_filled_groups():
    """
    Identifies and prints the nonabelian filled groups of order 2q^m for an odd prime q.

    The solution is based on the classification of finite filled groups. This script
    applies the constraints from the problem to the known list of filled groups.
    """
    print("This script identifies the nonabelian filled groups of order 2q^m,")
    print("where q is an odd prime and m is a natural number, based on the known classification theorem.\n")

    # Step 1: State the classification theorem for finite filled groups.
    # The source is P. Hegarty's paper "On filled groups" (2009).
    print("--- The Classification Theorem for Finite Filled Groups ---")
    print("A finite group G is filled if and only if it is one of the following:")
    print(" (i)   An elementary abelian 2-group (e.g., C2, C2 x C2).")
    print(" (ii)  A group of order p or p^2 for some odd prime p (e.g., C3, C5, C7, C9, C3 x C3).")
    print(" (iii) The cyclic group C4.")
    print(" (iv)  The dihedral groups D_6, D_8, D_10.")
    print(" (v)   The quaternion group Q_8.")
    print("----------------------------------------------------------\n")

    print("--- Analysis of the Conditions ---")
    print("We apply two conditions to this list:")
    print(" 1. The group must be nonabelian.")
    print(" 2. The group's order must be 2 * q^m for an odd prime q and m >= 1.\n")

    solutions = []

    # Check the list of filled groups
    
    # Cases (i), (ii), (iii) are all abelian groups, so we can discard them.
    # Now check the nonabelian groups from cases (iv) and (v).
    
    # Case D_6 (Dihedral group of order 6, also known as S_3)
    # It is nonabelian. Its order is 6.
    # We check if the order fits the form 2 * q^m.
    # 6 = 2 * 3. Here, 3 is an odd prime, so q=3 and m=1. This is a match.
    solutions.append({'group_name': 'D_6 (the Dihedral Group of order 6)', 'order': 6, 'q': 3, 'm': 1})

    # Case D_8 (Dihedral group of order 8)
    # It is nonabelian. Its order is 8.
    # We check if the order fits the form 2 * q^m.
    # 8 = 2 * 4. Here, 4 is not a power of an odd prime (4 = 2^2). No match.
    
    # Case D_10 (Dihedral group of order 10)
    # It is nonabelian. Its order is 10.
    # We check if the order fits the form 2 * q^m.
    # 10 = 2 * 5. Here, 5 is an odd prime, so q=5 and m=1. This is a match.
    solutions.append({'group_name': 'D_10 (the Dihedral Group of order 10)', 'order': 10, 'q': 5, 'm': 1})
    
    # Case Q_8 (Quaternion group of order 8)
    # It is nonabelian. Its order is 8.
    # We check if the order fits the form 2 * q^m.
    # 8 = 2 * 4. As with D_8, this does not match the required order form.
    
    print("--- Conclusion ---")
    if not solutions:
        print("There are no nonabelian filled groups whose order matches the form 2q^m.")
    else:
        print("The nonabelian filled groups of order 2q^m are:\n")
        for sol in solutions:
            print(f"Group: {sol['group_name']}")
            # Here we output each number in the final equation for the order.
            print(f"Equation of order: {sol['order']} = 2 * {sol['q']}^{sol['m']}\n")

if __name__ == '__main__':
    find_nonabelian_filled_groups()
