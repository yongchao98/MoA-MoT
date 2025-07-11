def find_smallest_n_for_reducible_space():
    """
    This script finds and explains the smallest non-negative integer n
    for which an n-point topological space can be reducible (not irreducible).

    A topological space X is reducible if it can be written as the union of
    a finite number of closed proper subsets. This is equivalent to being
    the union of two closed proper subsets, Z1 and Z2.
    """

    # We check for the smallest non-negative integer n.
    # n=0: The empty space. Its only closed set is itself (the empty set).
    #      It has no proper closed subsets, so it must be irreducible.
    # n=1: A one-point space X = {p}. The only possible topology is {emptyset, X}.
    #      The closed sets are also {emptyset, X}. The only proper closed subset
    #      is the empty set. Its union with itself is still the empty set, not X.
    #      So, any 1-point space is irreducible.
    # n=2: Let's consider a 2-point space X = {p1, p2}.
    #      We need to find a topology and two proper closed sets Z1, Z2
    #      such that Z1 U Z2 = X.
    n = 2

    print(f"The smallest non-negative integer for a reducible space is n = {n}.")
    print("\nHere is the construction for an n-point space where n = 2:")

    # Define the space and its points
    point_1 = "p1"
    point_2 = "p2"
    space_X = f"{{{point_1}, {point_2}}}"

    print(f"Let the topological space be X = {space_X}.")
    print("Consider the discrete topology on X, where every subset of X is both open and closed.")

    # Define the proper closed subsets
    closed_set_Z1 = f"{{{point_1}}}"
    closed_set_Z2 = f"{{{point_2}}}"

    print(f"Let Z1 = {closed_set_Z1}. This is a proper closed subset of X.")
    print(f"Let Z2 = {closed_set_Z2}. This is also a proper closed subset of X.")

    # Show that their union is X
    print("\nThe union of these two sets is:")
    print(f"Z1 U Z2 = {closed_set_Z1} U {closed_set_Z2} = {space_X}")
    print("Since Z1 U Z2 = X, the space X is reducible (not irreducible).")

    print(f"\nThus, the smallest such n is {n}.")
    
    # Return the final numerical answer
    return n

# Execute the function to print the explanation.
final_answer = find_smallest_n_for_reducible_space()
# The final answer is submitted in the special format below.
# The code above has already printed the detailed explanation and derivation.
