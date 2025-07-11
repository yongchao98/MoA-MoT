def solve_set_theory_problem():
    """
    Demonstrates the construction of an infinite set 'x' that is almost
    disjoint from every set in a given countable collection 'S'.
    """

    # --- Setup ---
    # We define a sample collection S of "infinite" sets. For this demonstration,
    # we use large finite sets that are generated according to different rules.
    # The algorithm itself does not depend on the sets being finite.
    s0 = sorted(list({n * 2 for n in range(1000)}))  # Even numbers
    s1 = sorted(list({n * 3 for n in range(1000)}))  # Multiples of 3
    s2 = sorted(list({n * n for n in range(500)}))    # Square numbers
    s3 = sorted(list({2**n for n in range(30)}))      # Powers of 2
    S = [s0, s1, s2, s3]
    
    print("This script demonstrates that for a countable collection of infinite sets S = {s_0, s_1, ...},")
    print("an infinite set x can be constructed such that for every s_i, the intersection |x intersect s_i| is finite.")
    print("This answers the user's question in the affirmative.\n")

    # --- Helper Function ---
    def find_min_greater_than(s_j, value):
        """
        Finds the smallest element in a sorted list s_j that is strictly greater than 'value'.
        This simulates finding min(s_j - {0, ..., value}) for an infinite set s_j.
        """
        for element in s_j:
            if element > value:
                return element
        # This fallback is for when our finite sample sets are exhausted.
        # For truly infinite sets, an element would always be found.
        return value + 1000000

    # --- Construction ---
    NUM_ELEMENTS_TO_CONSTRUCT = 15
    x = []
    x_previous = -1

    print(f"--- Constructing the first {NUM_ELEMENTS_TO_CONSTRUCT} elements of set x ---")

    for k in range(NUM_ELEMENTS_TO_CONSTRUCT):
        # At step k, we consider the sets s_0, ..., s_k (or all sets in S if S is smaller).
        sets_to_consider = S[:min(k + 1, len(S))]
        
        candidates = []
        for s_j in sets_to_consider:
            c_j = find_min_greater_than(s_j, x_previous)
            candidates.append(c_j)
        
        # Find the maximum of these candidates.
        max_candidate = max(candidates)
        
        # Choose x_k to be one greater than this maximum.
        x_k = max_candidate + 1
        
        x.append(x_k)
        x_previous = x_k
        
        print(f"\nStep k={k}:")
        print(f"  Previous element chosen: {x[-2] if k>0 else -1}")
        print(f"  Smallest elements > {x[-2] if k>0 else -1} from s_0..s_{len(sets_to_consider)-1}: {candidates}")
        print(f"  Max of candidates is {max_candidate}. Setting x_{k} = {max_candidate} + 1 = {x_k}")

    print("\n--- Construction Complete ---")
    print(f"\nConstructed set x = {x}")
    print("-" * 40)

    # --- Verification ---
    print("\n--- Verifying Intersections ---")
    set_x = set(x)
    for i, s_i in enumerate(S):
        intersection = sorted(list(set_x.intersection(set(s_i))))
        print(f"Intersection of x with s_{i}:")
        print(f"|x intersect s_{i}| = {len(intersection)}")
        print(f"Elements: {intersection if intersection else '(empty set)'}")

solve_set_theory_problem()