def solve():
    """
    Finds the smallest non-negative integer n for which an n-point
    topological space can be not irreducible (i.e., reducible).
    """

    print("--- Problem Definition ---")
    print("A topological space X is reducible if it can be written as a finite union of its proper closed subsets.")
    print("X = Z_1 U Z_2 U ... U Z_k, where each Z_i is a closed set and Z_i is a proper subset of X.")
    print("We are looking for the smallest 'n' (number of points) for which such a space exists.")
    print("\n--- Step-by-Step Analysis ---")

    print("\nCase n=0: The space is X = {}.")
    print("The only closed subset is X itself. There are no proper closed subsets.")
    print("Therefore, the 0-point space is always irreducible.")

    print("\nCase n=1: The space is X = {p}.")
    print("The only proper subset of X is the empty set {}. The union of any number of copies of {} is still {}, not X.")
    print("Therefore, any 1-point space is always irreducible.")

    print("\nCase n=2: Consider the space X = {1, 2}.")
    print("To show this space can be reducible, we need to find a topology where X is a union of its proper closed subsets.")
    print("Let's try to express X as the union of Z_1 = {1} and Z_2 = {2}.")
    print("This requires a topology where both {1} and {2} are closed sets.")
    print("Consider the discrete topology on X, where every subset is open. This implies every subset is also closed.")
    print("In this topology:")
    print(" - Z_1 = {1} is a proper closed subset of X.")
    print(" - Z_2 = {2} is also a proper closed subset of X.")
    print("\nNow, let's check their union, which forms the equation for reducibility:")
    
    z1_val = 1
    z2_val = 2
    
    print(f"X = {{{z1_val}}} U {{{z2_val}}}")
    print("The union is {1, 2}, which is the entire space X.")
    
    print("\n--- Conclusion ---")
    print("Since a 2-point space can be reducible, and spaces with 0 or 1 point cannot,")
    print("the smallest non-negative integer n is 2.")

    # The final answer
    answer = 2
    print(f"\n<<< {answer} >>>")

solve()