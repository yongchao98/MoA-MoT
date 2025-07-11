def solve_irreducible_space_problem():
    """
    This function explains and solves the problem of finding the smallest n
    for which a reducible n-point topological space exists.
    """
    print("Problem: Find the smallest non-negative integer n such that there exists an n-point topological space that is not irreducible.")
    
    print("\nStep 1: Definition of a Reducible Space")
    print("A topological space X is reducible if it can be written as a finite union of its proper closed subsets.")
    print("That is, X = Z_1 U Z_2 U ... U Z_k, where each Z_i is a closed set and Z_i is a proper subset of X.")

    print("\nStep 2: Analyzing cases for small n")
    
    print("\n--- Case n = 0 ---")
    print("Let X = {}. The only closed set is X itself. There are no proper closed subsets, so the 0-point space is irreducible.")

    print("\n--- Case n = 1 ---")
    print("Let X = {p}. The only proper closed subset is the empty set {}. A union of empty sets is still {}, not X. So, any 1-point space is irreducible.")

    print("\n--- Case n = 2 ---")
    print("Let X = {1, 2}. We need to find a topology where X is reducible.")
    print("Let's try to write X as a union of two proper closed subsets, Z_1 = {1} and Z_2 = {2}.")
    print("This requires a topology where both {1} and {2} are closed sets.")
    print("In the discrete topology, every subset is open, and therefore every subset is also closed.")
    print("With the discrete topology on X = {1, 2}:")
    print("  - Z_1 = {1} is a proper closed subset.")
    print("  - Z_2 = {2} is a proper closed subset.")
    print("We can write X as their union. Here is the equation showing the decomposition:")
    
    point1 = 1
    point2 = 2
    print(f"  X = {{{point1}}} U {{{point2}}}")

    print("\nStep 3: Conclusion")
    print("We have shown that for n=0 and n=1, all spaces are irreducible.")
    print("We found an example of a reducible space for n=2.")
    
    final_answer = 2
    print(f"\nTherefore, the smallest non-negative integer n is {final_answer}.")

solve_irreducible_space_problem()