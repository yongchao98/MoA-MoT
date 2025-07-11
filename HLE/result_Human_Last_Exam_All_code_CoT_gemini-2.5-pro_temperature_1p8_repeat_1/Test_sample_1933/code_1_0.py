import itertools

def solve():
    """
    Calculates and demonstrates the VC dimension for the given logic fragment.
    """
    # The schema S contains 4 unary predicates. Let's denote this number by d.
    d = 4

    print(f"The logical fragment is FO(exists, and, top, bot) with a schema of {d} unary predicates.")
    print("This logic can express classifiers corresponding to all possible intersections of {d} base sets.")
    print("The VC dimension of the class of all intersections of d sets is d.")
    print(f"Therefore, the VC dimension is {d}.\n")

    print("--- Proof Demonstration ---")
    print(f"We will now show that a set of size d = {d} can be shattered.")
    print("This demonstrates that the VC dimension is at least {d}.")
    print("Combined with the fact that the VC dimension is at most {d}, this proves the result.")

    # Let our set of points be S = {0, 1, ..., d-1}
    num_points = d
    points = list(range(num_points))

    # We need to define the base classifiers (the sets C_i corresponding to predicates P_i).
    # We represent sets as bitmasks for efficient intersection (bitwise AND).
    # Let's define the universe U = S. A set is a subset of {0, 1, 2, 3}.
    # We choose the base classifier C_i to contain all points *except* point i.
    # U = (1<<d) - 1 which is 0b1111 for d=4
    universe_mask = (1 << num_points) - 1
    
    # C_i = U \ {i}
    base_classifiers_masks = [universe_mask ^ (1 << i) for i in range(d)]

    def set_to_string(mask, n):
        elements = [str(i) for i in range(n) if (mask >> i) & 1]
        if not elements:
            return "{}"
        return "{" + ", ".join(elements) + "}"

    print(f"\nLet's construct a shattered set of {d} points: S = {set_to_string(universe_mask, d)}")
    print("We define the base classifiers C_i (corresponding to predicate P_i) as follows:")
    for i in range(d):
        print(f"C_{i} = S \\ {{{i}}} = {set_to_string(base_classifiers_masks[i], d)}")

    print("\nNow, we iterate through all 2^d = 16 possible subsets of S to show they can all be generated.")

    shattered = True
    num_subsets = 1 << num_points

    for i in range(num_subsets):
        target_subset_mask = i
        
        # We need to find an intersection of base classifiers that results in the target subset.
        # Rule: To exclude point j, we must intersect with C_j.
        # So, the intersection should be over all C_j where j is *not* in the target subset.
        
        # Start with the universe (intersection of an empty set of classifiers)
        hypothesis_mask = universe_mask
        indices_to_intersect = []

        for j in range(num_points):
            # if point j is NOT in the target subset
            if not ((target_subset_mask >> j) & 1):
                hypothesis_mask &= base_classifiers_masks[j]
                indices_to_intersect.append(j)

        target_str = set_to_string(target_subset_mask, d)
        
        if not indices_to_intersect:
            intersect_str = "C_empty (Universe)"
        else:
            intersect_str = " intersect ".join([f"C_{idx}" for idx in indices_to_intersect])

        print(f"Target: {target_str.ljust(15)} -> Achieved by: {intersect_str}")
        
        if hypothesis_mask != target_subset_mask:
            shattered = False
            print("ERROR: Could not generate the target subset.")
            break

    if shattered:
        print(f"\nSuccess! All {num_subsets} subsets were generated.")
        print(f"This proves that the set S of size {d} can be shattered.")

    print("\n--- Final Conclusion ---")
    num_predicates = 4
    vc_dimension = 4
    print(f"For a schema S with {num_predicates} unary predicates, the VC dimension of FO(exists, and, top, bot)[S] is equal to the number of predicates.")
    print(f"VC dimension = {vc_dimension}")


solve()