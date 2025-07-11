import sys

def solve():
    """
    This function finds the smallest non-negative integer n such that there exists
    an n-point topological space that is not irreducible.
    """

    print("Analyzing the problem step-by-step:")
    print("-" * 30)

    # Step 1: Analyze for n=0
    print("Case n=0: The space is X = {}.")
    print("The only subset (and only closed subset) of X is {} itself.")
    print("A proper subset must be different from the space X.")
    print("Therefore, X has no proper closed subsets and cannot be written as their union.")
    print("Result: Any 0-point space is irreducible.\n")

    # Step 2: Analyze for n=1
    print("Case n=1: The space is X = {p} (a single point).")
    print("The proper subsets of X are just {}.")
    print("The empty set {} is always a closed set.")
    print("Can we form X by unioning proper closed subsets? No, the union would just be {}.")
    print("Result: Any 1-point space is irreducible.\n")

    # Step 3: Analyze for n=2
    print("Case n=2: The space is X = {a, b}.")
    print("To show this space is NOT irreducible, we need to write X as a union of proper closed subsets.")
    print("The proper subsets of X are {}, {a}, and {b}.")
    print("The only way to form X by union is: {a} U {b} = X.")
    print("This requires finding a topology on X where both {a} and {b} are closed sets.")
    print("A set is closed if its complement is open.")
    print("  - For {a} to be closed, its complement X \\ {a} = {b} must be open.")
    print("  - For {b} to be closed, its complement X \\ {b} = {a} must be open.")
    print("The discrete topology T = { {}, {a}, {b}, {a, b} } satisfies these conditions.")
    print("With this topology, {a} and {b} are indeed proper closed subsets.")
    print("Let's demonstrate the union:")

    # Define the sets for the equation
    z1 = {'a'}
    z2 = {'b'}
    x_space = z1.union(z2)
    
    # Python's set printing uses repr, which gives {'a'} etc.
    # We want a more readable output for sets like {a}.
    set_to_str = lambda s: '{' + ', '.join(sorted(list(s))) + '}' if s else '{}'

    # The required final equation outputting each number (or element)
    print(f"\nFinal Equation: {set_to_str(z1)} U {set_to_str(z2)} = {set_to_str(x_space)}\n")
    
    print("This confirms that a 2-point space can be non-irreducible.")
    print("Since spaces with n=0 and n=1 must be irreducible, the smallest such integer is 2.")

    # Final Answer
    final_answer = 2
    sys.stdout.write(f"\n<<<{final_answer}>>>\n")

solve()