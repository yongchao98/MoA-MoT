import itertools

def get_powerset(s):
    """Returns the powerset of a given set."""
    # The chain.from_iterable is a fast way to flatten a list of lists.
    # combinations(s, r) generates all subsets of size r.
    # We iterate r from 0 to len(s) to get all possible subset sizes.
    s_list = list(s)
    return set(itertools.chain.from_iterable(itertools.combinations(s_list, r) for r in range(len(s_list) + 1)))

def solve():
    """
    Calculates and demonstrates the VC dimension of the specified logic fragment.
    """
    # Number of unary predicates in the schema S
    num_predicates = 4
    
    # The VC dimension of the class of intersections of k sets is k.
    vc_dimension = num_predicates

    print(f"The schema S has {num_predicates} unary predicates.")
    print("The concept class consists of all possible intersections of the sets defined by these predicates.")
    print("This is known as the class of monotone monomials over 4 variables.")
    print("The VC dimension of this class is equal to the number of predicates.\n")

    # --- Demonstration that VC dimension is at least 4 ---
    print(f"To demonstrate that the VC dimension is at least {num_predicates}, we will show that a set of {num_predicates} points can be shattered.")
    
    # A set of 4 points to be shattered
    X = frozenset(range(num_predicates))
    print(f"Let's define a set of {num_predicates} points, X = {set(X)}")
    
    # We define 4 base hypotheses H_i.
    # The hypothesis H_i is the set of all points in X except for point i.
    # This construction is a standard way to prove the lower bound on the VC-dim for monomials.
    base_hypotheses = [X - {i} for i in X]
    
    print("\nWe define the 4 base concepts (interpretations of P_i(x)) as follows:")
    for i, h in enumerate(base_hypotheses):
        print(f"H_{i+1} = {set(h)}")

    # Generate the powerset of X, which are all the subsets we need to form.
    powerset_X = get_powerset(X)
    
    print(f"\nNow, we will show that we can generate all {2**num_predicates} subsets of X by taking intersections of these base hypotheses.")
    
    generated_subsets = set()
    
    # Iterate through all possible subsets of the base hypotheses to form all possible intersections
    for i in range(len(base_hypotheses) + 1):
        for hs_to_intersect_indices in itertools.combinations(range(len(base_hypotheses)), i):
            
            # Start with the universe (X) for the intersection
            current_intersection = set(X)
            
            # The intersection of an empty set of hypotheses is the universe
            if not hs_to_intersect_indices:
                concept = set(X)
            else:
                # Intersect the selected hypotheses
                hypotheses_to_intersect = [base_hypotheses[j] for j in hs_to_intersect_indices]
                concept = set.intersection(*map(set, hypotheses_to_intersect))
                
            generated_subsets.add(frozenset(concept))

    # Verify that all subsets of X were generated.
    if generated_subsets == powerset_X:
        print(f"\nSuccess! All {len(powerset_X)} subsets of X were successfully generated.")
        print("This shatters the set X and confirms that the VC dimension is at least 4.")
    else:
        print("\nFailure! Not all subsets could be generated.")
        print(f"Missing subsets: {powerset_X - generated_subsets}")

    print("\nSince the VC dimension of k monomials is exactly k, and we have shown it's at least 4, we can conclude:")
    # Final output as requested
    equation_parts = ["VC", "dimension", "=", str(vc_dimension)]
    for part in equation_parts:
        print(part, end=" ")
    print()


solve()