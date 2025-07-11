def find_smallest_n():
    """
    This script explains the reasoning to find the smallest non-negative integer n
    such that an n-point topological space can be reducible (not irreducible).
    """

    print("--- Problem Definition ---")
    print("A topological space X is 'not irreducible' (or reducible) if it can be written as a union of finitely many proper closed subsets.")
    print("A 'proper' subset is a subset that is not equal to the whole space.")
    print("We are looking for the smallest non-negative integer n for which such a space exists.")

    print("\n--- Analysis of n=0 and n=1 ---")
    print("For n=0 (the empty set), there are no proper subsets at all, so it must be irreducible.")
    print("For n=1 (a single point space), the only proper subset is the empty set, which must be closed. A union of empty sets is still the empty set, not the whole space. So, it must be irreducible.")
    
    print("\n--- Analysis of n=2 ---")
    print("Let's consider a 2-point space, X = {0, 1}.")
    print("To show it can be reducible, we need to find a topology where X is a union of its proper closed subsets.")
    print("The only way to form X from its proper subsets is as the union {0} U {1}.")
    print("This requires that both the set {0} and the set {1} are closed sets in our topology.")
    
    print("\nConsider the discrete topology on X, where all subsets are open.")
    print("Open sets: {∅, {0}, {1}, {0, 1}}")
    print("The closed sets are their complements: {{0, 1}, {1}, {0}, ∅}")

    print("\nWe can now demonstrate that X is reducible:")
    print("Let the first proper closed set be Z1 = {0}.")
    print("Let the second proper closed set be Z2 = {1}.")
    print("Their union is Z1 U Z2 = X.")
    
    print("\nThe equation showing the decomposition is: {0, 1} = {0} U {1}")
    
    # Printing the numbers in the final equation as requested
    z1_elements = [0]
    z2_elements = [1]
    
    for number in z1_elements:
        print("The number in the first proper closed set is:", number)
    
    for number in z2_elements:
        print("The number in the second proper closed set is:", number)
        
    print("\n--- Conclusion ---")
    print("Since n=0 and n=1 spaces are always irreducible, and we have constructed a reducible space for n=2, the smallest such integer is 2.")


if __name__ == "__main__":
    find_smallest_n()