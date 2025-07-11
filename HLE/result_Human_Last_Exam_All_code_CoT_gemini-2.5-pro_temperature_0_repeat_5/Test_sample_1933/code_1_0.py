import itertools

def get_shattered_set(k):
    """
    Generates a set of k points that can be shattered by monotone conjunctions.
    Each point is a tuple of length k with k-1 ones and one zero.
    """
    points = []
    for i in range(k):
        point = [1] * k
        point[i] = 0
        points.append(tuple(point))
    return points

def check_conjunction(point, predicate_indices):
    """
    Checks if a point satisfies a monotone conjunction.
    The conjunction is defined by the indices of the predicates.
    A point (x_1, ..., x_k) satisfies P_i1 AND P_i2 ...
    if x_i1 = 1, x_i2 = 1, ...
    """
    for i in predicate_indices:
        # In programming, indices are 0-based, so we subtract 1
        if point[i - 1] == 0:
            return False
    return True

def find_shattering_formula(subset, all_points, k):
    """
    Finds the monotone conjunction that selects the given subset of points.
    """
    if not subset:
        # To select the empty set, we need a conjunction that is false for all points.
        # The conjunction of all k predicates works.
        return list(range(1, k + 1))

    # The set of indices of points in the target subset
    subset_indices = {all_points.index(p) for p in subset}

    # The required predicate indices are the complement of the point indices
    # in the subset.
    # e.g., to select points {p_0, p_2}, we need predicates P_1, P_3 (if k=4)
    predicate_indices = [i + 1 for i in range(k) if i not in subset_indices]
    return predicate_indices

def main():
    """
    Demonstrates the shattering of k points for k=4.
    """
    k = 4
    print(f"Let S be a schema with k={k} unary predicates: P_1, P_2, P_3, P_4.")
    print("We show that the VC dimension is at least 4 by shattering a set of 4 points.")
    
    # 1. Define the set of k points to be shattered
    points_to_shatter = get_shattered_set(k)
    print(f"\nConsider the following set of {k} points (or elements):")
    for i, p in enumerate(points_to_shatter):
        print(f"  p{i+1} = {p}  (Represents a type where P_j is true if the j-th component is 1)")

    # 2. Iterate through all 2^k subsets of these points
    print("\nWe will now show that any subset of these points can be selected by a monotone conjunction.")
    
    num_subsets = 0
    for i in range(len(points_to_shatter) + 1):
        for subset in itertools.combinations(points_to_shatter, i):
            num_subsets += 1
            subset = set(subset)
            
            # 3. For each subset, find the conjunction that selects it
            predicate_indices = find_shattering_formula(subset, points_to_shatter, k)
            
            # 4. Format and print the result
            subset_str = "{" + ", ".join(sorted([f"p{points_to_shatter.index(p)+1}" for p in subset])) + "}"
            if not predicate_indices:
                formula_str = "TOP (the empty conjunction)"
            else:
                formula_str = " AND ".join([f"P_{j}(x)" for j in predicate_indices])
            
            print(f"\n- Subset {subset_str}:")
            print(f"  Is selected by the formula: {formula_str}")

            # 5. Verification step
            selected_points = {p for p in points_to_shatter if check_conjunction(p, predicate_indices)}
            if selected_points == subset:
                print("  Verification: Correct.")
            else:
                print("  Verification: Incorrect.")

    print(f"\nSuccessfully demonstrated that all {num_subsets} subsets can be generated.")
    print("This shows that the VC dimension is at least 4.")
    print("Since there are only 2^4 = 16 monotone conjunctions, we cannot shatter more than 4 points.")
    print("Thus, the VC dimension is exactly 4.")

if __name__ == "__main__":
    main()