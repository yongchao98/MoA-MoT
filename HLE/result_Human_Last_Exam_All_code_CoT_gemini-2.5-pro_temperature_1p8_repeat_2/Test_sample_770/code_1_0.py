import sympy

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.

    This rank is equal to the number of non-identity conjugacy classes of the group A_5.
    """
    # The orientation-preserving icosahedral group is isomorphic to the Alternating Group on 5 elements, A_5.
    G = sympy.combinatorics.group.AlternatingGroup(5)

    # The number of conjugacy classes in a group can be computed directly.
    # For A_5, these correspond to cycle types: id, (2,2), (3), and two classes of (5)-cycles.
    num_total_classes = G.num_conjugacy_classes

    # The identity element always forms its own conjugacy class.
    # The number we are interested in is the number of non-identity classes.
    num_identity_classes = 1

    # The rank of H^2_c(Y, Q) is the number of non-identity conjugacy classes.
    rank = num_total_classes - num_identity_classes

    print("Step 1: Identify the group G.")
    print("The orientation-preserving icosahedral group is G = A_5.")
    print("-" * 20)
    print("Step 2: Relate cohomology to the group structure via the McKay correspondence.")
    print("The rank of H^2_c(Y, Q) is the number of non-identity conjugacy classes of G.")
    print("-" * 20)
    print("Step 3: Calculate the number of classes.")
    print(f"The total number of conjugacy classes in A_5 is: {num_total_classes}")
    print(f"The number of identity classes in A_5 is: {num_identity_classes}")
    print(f"The rank of H^2_c(Y, Q) is the difference:")
    print(f"{num_total_classes} - {num_identity_classes} = {rank}")

solve_cohomology_rank()
<<<4>>>