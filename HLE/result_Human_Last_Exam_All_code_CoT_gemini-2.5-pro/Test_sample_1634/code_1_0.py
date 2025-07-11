import sys

def solve():
    """
    This function determines the smallest non-negative integer n such that there exists
    an n-point topological space that is not irreducible.
    It prints the step-by-step reasoning.
    """
    print("Goal: Find the smallest non-negative integer n such that an n-point topological space exists that is NOT irreducible.")
    print("\nDefinition: A topological space X is reducible (not irreducible) if it can be written as the union of two proper closed subsets.")
    print("Let's say X = Z1 U Z2, where Z1 and Z2 are closed and Z1 != X, Z2 != X.")

    print("\n--- Analysis for n = 0 ---")
    print("Let X be a 0-point space, so X = {}.")
    print("The only topology is T = {{}}.")
    print("The closed sets are the complements of open sets. The only closed set is X itself.")
    print("There are no *proper* closed subsets. So X cannot be written as a union of proper closed subsets.")
    print("Result: Any 0-point space is irreducible.")

    print("\n--- Analysis for n = 1 ---")
    print("Let X be a 1-point space, so X = {p}.")
    print("The proper subsets of X are just {}.")
    print("In any topology, {} is a closed set. This is the only possible proper closed set.")
    print("A union of proper closed sets can only be {}, which is not equal to X.")
    print("Result: Any 1-point space is irreducible.")

    print("\n--- Analysis for n = 2 ---")
    print("Let X be a 2-point space, so X = {a, b}.")
    print("We want to find a topology on X such that X is reducible.")
    print("This requires X = Z1 U Z2, for proper closed subsets Z1 and Z2.")
    print("To get the union to be X, we must have Z1 = {a} and Z2 = {b} (or vice-versa).")
    print("So, we need a topology where both {a} and {b} are closed sets.")
    print("Consider the discrete topology on X, where every subset is open.")
    print("   Open sets: {}, {a}, {b}, {a, b}")
    print("   Closed sets (complements): {a, b}, {b}, {a}, {}")
    print("Let's check if this space is reducible:")
    print("   - Let Z1 = {a}. This is a proper closed subset.")
    print("   - Let Z2 = {b}. This is also a proper closed subset.")
    print("Now, let's look at the union, which represents our equation:")
    print("   Equation: Z1 U Z2 = X")
    # The prompt requests to "output each number in the final equation"
    # Here, we represent the elements of the sets.
    print("   Substituting the sets: {a} U {b} = {a, b}")
    print("The union is indeed X. Therefore, this 2-point space is reducible.")
    
    print("\n--- Conclusion ---")
    final_answer = 2
    print(f"We found that for n=0 and n=1, all spaces are irreducible.")
    print(f"For n={final_answer}, we constructed an example of a reducible space.")
    print(f"The smallest non-negative integer n is {final_answer}.")
    
    # This part is for the platform to capture the final numerical answer.
    # It is not intended for the user to see in a normal run.
    # sys.stdout.write(f"\n<<<{final_answer}>>>")


solve()