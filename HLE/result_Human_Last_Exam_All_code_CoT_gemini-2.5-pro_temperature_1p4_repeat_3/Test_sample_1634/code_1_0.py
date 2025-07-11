def find_smallest_n():
    """
    This script determines the smallest non-negative integer n for which
    an n-point topological space that is not irreducible (i.e., is reducible) can exist.
    """
    print("We are looking for the smallest integer n >= 0 such that an n-point space can be 'reducible'.")
    print("A space X is reducible if it can be written as a finite union of its proper closed subsets.")
    print("-" * 40)

    # Case n = 0
    print("Case n = 0: The space is X = {}.")
    print("The only subset is {}, which is not a 'proper' subset of X.")
    print("Thus, the 0-point space is always irreducible.")
    print("-" * 40)

    # Case n = 1
    print("Case n = 1: The space is X = {p}.")
    print("The only proper subset of X is the empty set, {}.")
    print("In any topology, the only possible proper closed subset is {}.")
    print("The union of empty sets cannot form X. Thus, the 1-point space is always irreducible.")
    print("-" * 40)

    # Case n = 2
    print("Case n = 2: The space is X = {0, 1}.")
    print("Let's see if we can define a topology that makes this space reducible.")
    print("Consider the discrete topology, where every subset is closed.")
    print("The proper closed subsets are {}, {0}, and {1}.")
    print("We can choose the proper closed subsets Z1 = {0} and Z2 = {1}.")
    print("Let's check their union: Z1 U Z2 = {0} U {1} = {0, 1}.")
    print("This union is equal to the entire space X.")

    print("\nWe have successfully represented X as a union of its proper closed subsets.")
    print("The final equation is:")

    # Printing each number in the final equation as requested.
    point1 = 0
    point2 = 1
    print(f"{{{point1}, {point2}}} = {{{point1}}} U {{{point2}}}")

    smallest_n = 2
    print(f"\nConclusion: A reducible space can exist for n = 2.")
    print(f"Since spaces with n=0 and n=1 are always irreducible, the smallest such non-negative integer is {smallest_n}.")

find_smallest_n()
<<<2>>>