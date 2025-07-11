def solve():
    """
    Finds and explains the smallest n for a non-irreducible n-point space.
    """
    print("This script determines the smallest non-negative integer n such that an n-point topological space that is not irreducible exists.")
    print("-" * 70)

    print("\n[Step 1: Definition of a Non-Irreducible Space]")
    print("A topological space X is defined as NOT irreducible (or reducible) if it can be written as a finite union of its proper closed subsets.")
    print("That is, X = Z1 U Z2 U ... U Zk, where each Zi is a closed set and Zi is a proper subset of X (Zi != X).")

    print("\n[Step 2: Analyzing n = 0]")
    print("Let X be a 0-point space, i.e., X = ∅.")
    print("The only closed set in this space is X itself (the empty set). There are no proper closed subsets.")
    print("Therefore, X cannot be written as a union of proper closed subsets. The 0-point space is irreducible.")

    print("\n[Step 3: Analyzing n = 1]")
    print("Let X be a 1-point space, i.e., X = {p}.")
    print("The closed sets are {p} and ∅. The only proper closed subset is ∅.")
    print("Any union of this subset is still ∅, which is not equal to X.")
    print("Therefore, any 1-point space is irreducible.")

    print("\n[Step 4: Analyzing n = 2]")
    print("Let X be a 2-point space, i.e., X = {p, q}.")
    print("We can define the discrete topology on X, where every subset is open.")
    print("In the discrete topology, every subset is also closed.")
    print("The proper subsets of X are ∅, {p}, and {q}. All of these are closed.")
    print("Let's choose two proper closed subsets: Z1 = {p} and Z2 = {q}.")
    print("Their union is Z1 U Z2 = {p} U {q} = {p, q}, which is equal to X.")
    print("\nHere is the decomposition showing X is not irreducible:")
    
    # The "equation" part requested by the user prompt
    space = "{p, q}"
    subset1 = "{p}"
    subset2 = "{q}"
    print(f"Space X = {space}")
    print(f"Proper closed subset Z1 = {subset1}")
    print(f"Proper closed subset Z2 = {subset2}")
    print(f"Union: {space} = {subset1} U {subset2}")

    print("\nSince we found an example for n=2, and n=0 and n=1 did not work, the smallest such integer is 2.")
    print("-" * 70)

solve()
<<<2>>>