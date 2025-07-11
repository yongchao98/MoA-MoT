import sys

def solve_irreducible_space_problem():
    """
    This function explains the step-by-step reasoning to find the smallest
    non-negative integer n for which an n-point topological space can be reducible.
    """
    
    print("Goal: Find the smallest non-negative integer n such that there exists an n-point topological space that is not irreducible.")
    
    print("\n--- Definition ---")
    print("A topological space X is irreducible if it is not a union of a finite number of its proper closed subsets.")
    print("A space is reducible (i.e., not irreducible) if it can be written as X = Z_1 U Z_2 U ... U Z_k, where each Z_i is a closed subset and Z_i is a proper subset of X (Z_i != X).")

    print("\n--- Case n = 0 ---")
    print("Let X be a 0-point space. Then X is the empty set, X = {}.")
    print("The only topology is T = {{}}. The only closed set is X itself.")
    print("There are no *proper* closed subsets. Therefore, X cannot be written as a union of proper closed subsets.")
    print("Conclusion: The 0-point space is irreducible.")

    print("\n--- Case n = 1 ---")
    print("Let X be a 1-point space, X = {p}.")
    print("Regardless of the topology, the set of closed sets must contain X and the empty set {}.")
    print("The only possible proper closed subset is the empty set {}.")
    print("Any finite union of the empty set is still the empty set, which is not equal to X.")
    print("Conclusion: Any 1-point space is irreducible.")

    print("\n--- Case n = 2 ---")
    print("Let X be a 2-point space. We will use the points {0, 1}, so X = {0, 1}.")
    print("We need to check if there *exists* a topology on X that makes it reducible.")
    print("Let's consider the discrete topology, where every subset is an open set.")
    print("The open sets are: {}, {0}, {1}, {0, 1}.")
    print("In any topology, a set is closed if its complement is open. In the discrete topology, this means every subset is also a closed set.")
    print("The closed sets are: {}, {0}, {1}, {0, 1}.")
    print("The proper closed subsets of X are those not equal to X. They are: {}, {0}, {1}.")
    print("Now, let's see if we can write X as a union of its proper closed subsets.")
    print("Let's choose two proper closed subsets: Z_1 = {0} and Z_2 = {1}.")
    print("The union of these two sets is Z_1 U Z_2 = {0} U {1} = {0, 1}.")
    print("\nThis gives us the final equation showing reducibility:")
    print("X = Z_1 U Z_2")
    # Using python sets to demonstrate.
    point_0 = 0
    point_1 = 1
    Z_1 = {point_0}
    Z_2 = {point_1}
    X = Z_1.union(Z_2)
    print(f"{X} = {{{point_0}}} U {{{point_1}}}")
    print("\nSince X is a union of two of its proper closed subsets, this space is reducible (not irreducible).")

    print("\n--- Final Conclusion ---")
    print("We have shown that for n=0 and n=1, any space is irreducible.")
    print("For n=2, we have constructed a space that is reducible.")
    print("Therefore, the smallest non-negative integer n such that there exists an n-point topological space that is not irreducible is 2.")

if __name__ == "__main__":
    solve_irreducible_space_problem()