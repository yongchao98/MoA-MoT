def find_smallest_n_for_reducible_space():
    """
    This script finds and explains the smallest non-negative integer n for which
    an n-point topological space can be not irreducible (i.e., reducible).

    A topological space X is reducible if it can be written as a finite
    union of proper closed subsets (where each subset Z_i is not equal to X).
    The condition is: X = Z_1 U Z_2 U ... U Z_k.
    """
    print("This script will determine the smallest non-negative integer n such that an n-point space can be reducible.")
    print("-" * 80)
    
    # Step 1: Analyze cases where n < 2
    print("Step 1: Analyze small values of n.")
    
    print("\n[Case n=0]")
    print("The space is X = {}, the empty set.")
    print("The only closed subset of X is X itself. There are no *proper* closed subsets.")
    print("Conclusion: The 0-point space is always irreducible.")

    print("\n[Case n=1]")
    print("The space is X = {p1}, a single point.")
    print("The only proper closed subset in any topology on X is the empty set, {}.")
    print("A union of empty sets is {}, which cannot equal X.")
    print("Conclusion: Any 1-point space is always irreducible.")

    # Step 2: Analyze the case n = 2
    print("\n[Case n=2]")
    print("Let the space be X = {1, 2}.")
    print("To show it can be reducible, we must find a topology where X is a union of proper closed subsets.")
    print("Let's consider the discrete topology on X. In this topology, every subset is open, and therefore every subset is also closed.")
    print("The set of all closed subsets is {{}, {1}, {2}, {1, 2}}.")
    print("The proper closed subsets are {}, {1}, and {2}.")
    
    # Define the sets for the equation
    Z1 = {1}
    Z2 = {2}
    X = Z1.union(Z2)
    
    print(f"\nWe can choose two proper closed subsets: Z1 = {Z1} and Z2 = {Z2}.")
    print("Their union is Z1 U Z2 = {1} U {2} = {1, 2}, which is the entire space X.")
    print("This decomposition shows that the space is reducible.")
    
    # Print the final equation with each number shown
    print("\nThe equation demonstrating this is:")
    print(f"  X      =    Z1    U    Z2")
    print(f" {X}   =   {Z1}   U   {Z2}")
    
    print("-" * 80)

    # Final Conclusion
    n = 2
    print("Since spaces for n=0 and n=1 are always irreducible, and we have constructed a reducible space for n=2,")
    print(f"the smallest such non-negative integer n is {n}.")

# Execute the reasoning
find_smallest_n_for_reducible_space()