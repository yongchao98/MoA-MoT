import sympy

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y
    of the quotient singularity C^3/G, where G is the icosahedral group.

    This rank is equal to the number of non-trivial conjugacy classes of G.
    The icosahedral group is isomorphic to the alternating group A_5.
    """
    # The icosahedral group is isomorphic to the alternating group on 5 elements, A_5.
    n = 5
    group_name = f"A_{n}"
    G = sympy.combinatorics.group.AlternatingGroup(n)

    # The number of conjugacy classes can be found directly.
    # For A_5, the classes are:
    # 1. The identity (1 element)
    # 2. Products of two disjoint transpositions, e.g., (1,2)(3,4) (15 elements)
    # 3. 3-cycles, e.g., (1,2,3) (20 elements)
    # 4. 5-cycles of one type, e.g., (1,2,3,4,5) (12 elements)
    # 5. 5-cycles of another type, e.g., (1,3,5,2,4) (12 elements)
    conjugacy_classes_list = G.conjugacy_classes()
    total_classes = len(conjugacy_classes_list)
    trivial_classes = 1 # The class of the identity element

    # The rank of H^2_c(Y, Q) is the number of non-trivial conjugacy classes.
    rank = total_classes - trivial_classes

    print(f"The group G is the icosahedral group, isomorphic to {group_name}.")
    print(f"The number of conjugacy classes in {group_name} is: {total_classes}")
    print(f"The number of trivial conjugacy classes is: {trivial_classes}")
    print("-" * 30)
    print("The rank of H^2_c(Y, Q) is the number of non-trivial conjugacy classes.")
    print(f"Rank = (Total Classes) - (Trivial Classes)")
    # The final print statement outputs each number in the equation, as requested.
    print(f"Rank of H^2_c(Y, Q) = {total_classes} - {trivial_classes} = {rank}")

solve_cohomology_rank()
<<<4>>>