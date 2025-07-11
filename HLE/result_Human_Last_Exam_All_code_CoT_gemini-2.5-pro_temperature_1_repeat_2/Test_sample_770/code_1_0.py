from sympy.combinatorics import AlternatingGroup

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/G,
    where G is the icosahedral group.

    The rank is equal to the number of non-trivial conjugacy classes of G.
    The icosahedral group G is isomorphic to the alternating group A_5.
    """
    # Step 1: Define the group G as the Alternating Group on 5 elements.
    G = AlternatingGroup(5)
    group_name = "A_5 (the alternating group on 5 elements)"
    print(f"The group G is the icosahedral group, which is isomorphic to {group_name}.")
    print("-" * 30)

    # Step 2: Get the conjugacy classes of the group.
    # The `conjugacy_classes` method returns a list of classes, where each class is a set of its elements.
    conjugacy_classes = G.conjugacy_classes()

    # Step 3: Count the total number of conjugacy classes.
    num_total_classes = len(conjugacy_classes)
    print(f"The total number of conjugacy classes in {group_name} is: {num_total_classes}")

    # Step 4: The number of non-trivial conjugacy classes is the total number minus one
    # (for the class of the identity element).
    # This number corresponds to the rank we are looking for.
    identity_classes = 1
    rank = num_total_classes - identity_classes
    
    print(f"The number of non-trivial conjugacy classes is the total number of classes minus the number of trivial classes (which is 1).")
    print(f"The final equation is: {num_total_classes} - {identity_classes} = {rank}")
    print("-" * 30)
    print(f"According to the McKay correspondence for dimension 3, this number is the rank of H^2_c(Y, Q).")
    print(f"The rank of H^2_c(Y, Q) is: {rank}")

if __name__ == '__main__':
    solve_cohomology_rank()