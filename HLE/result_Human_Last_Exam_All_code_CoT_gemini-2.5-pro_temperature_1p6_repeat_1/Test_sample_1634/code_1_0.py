def find_smallest_n_for_reducible_space():
    """
    This function explains the step-by-step reasoning to find the smallest
    non-negative integer n such that there exists an n-point topological space
    that is not irreducible (i.e., is reducible).
    """

    print("--- Problem Definition ---")
    print("A topological space X is irreducible if it cannot be written as a finite union of proper closed subsets.")
    print("A space is NOT irreducible (i.e., is reducible) if it can be written as X = Z1 U Z2,")
    print("where Z1 and Z2 are proper closed subsets of X (Z1 != X and Z2 != X).")
    print("We are looking for the smallest non-negative integer 'n' for which an n-point space can be reducible.")
    print("\n--- Analysis ---")

    # Case n = 0
    print("\nStep 1: Checking n = 0")
    print("Let X be a 0-point space, so X = {}.")
    print("The only closed set in this space is X itself. There are no *proper* closed subsets.")
    print("Therefore, X cannot be written as a union of proper closed subsets.")
    print("Result: The 0-point space is irreducible.")

    # Case n = 1
    print("\nStep 2: Checking n = 1")
    print("Let X be a 1-point space, so X = {p}.")
    print("The only proper subset of X is the empty set {}.")
    print("In any topology on X, the empty set is always a closed set.")
    print("To be reducible, X would need to be a union of proper closed subsets. However, the union of empty sets is still the empty set, not X.")
    print("Result: Any 1-point space is irreducible.")

    # Case n = 2
    print("\nStep 3: Checking n = 2")
    print("Let X be a 2-point space, so X = {p1, p2}.")
    print("To show this space can be reducible, we need to find a topology on X that allows X = Z1 U Z2, where Z1 and Z2 are proper closed subsets.")
    print("Let's try to set Z1 = {p1} and Z2 = {p2}.")
    print("For Z1 and Z2 to be closed sets, their complements must be open. This is true in the discrete topology, where every subset is open (and also closed).")
    print("In the 2-point discrete space, the sets Z1 = {p1} and Z2 = {p2} are proper closed subsets.")
    print("Now we check if their union forms X, based on the equation Z1 U Z2 = X:")
    point_1 = "p1"
    point_2 = "p2"
    z1 = {point_1}
    z2 = {point_2}
    x_set = z1.union(z2)
    
    # Printing the sets in the final equation explicitly
    print(f"The set Z1 is {{{next(iter(z1))}}}.")
    print(f"The set Z2 is {{{next(iter(z2))}}}.")
    print(f"Their union is {{{', '.join(sorted(list(x_set)))}}}.")
    print(f"So the equation is: {{{next(iter(z1))}}} U {{{next(iter(z2))}}} = {{{', '.join(sorted(list(x_set)))}}}")
    print("This shows that X is the union of two of its proper closed subsets.")
    print("Result: A 2-point space with the discrete topology is NOT irreducible.")

    # Conclusion
    print("\n--- Conclusion ---")
    print("We have shown:")
    print("- For n=0, any space is irreducible.")
    print("- For n=1, any space is irreducible.")
    print("- For n=2, a reducible space exists.")
    print("\nTherefore, the smallest non-negative integer n is 2.")

find_smallest_n_for_reducible_space()