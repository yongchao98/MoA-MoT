import itertools

def demonstrate_shattering():
    """
    This function demonstrates that the VC dimension for a schema with 4 unary
    predicates is at least 4 by showing a set of 4 points that can be shattered.
    """
    k = 4
    print(f"Let the number of unary predicates be k = {k}.")
    print("We will show that a set of k points can be shattered.")
    print("-" * 20)

    points = [f"a{i+1}" for i in range(k)]
    predicates = [f"P{j+1}" for j in range(k)]

    # We define a model (the interpretation of predicates on the points).
    # P_j(a_i) is true if and only if i != j.
    # model[i][j] represents the truth value of P_{j+1}(a_{i+1}).
    model = [[(i != j) for j in range(k)] for i in range(k)]

    print("Consider a set of 4 points: " + ", ".join(points))
    print("We define the 4 predicates over these points as follows:")
    print("P_j(a_i) is TRUE if i != j, and FALSE if i == j.")
    print("Predicate truth table:")
    header = "      " + " ".join([f"{p:<5}" for p in predicates])
    print(header)
    print("      " + "----- " * k)
    for i in range(k):
        row_str = f"{points[i]:<5} |" + "".join([f" {str(model[i][j]):<5}" for j in range(k)])
        print(row_str)
    print("-" * 20)


    shattered = True
    num_subsets = 2**k
    print(f"Now, let's check all {num_subsets} possible subsets (dichotomies) of the points.")

    for i in range(num_subsets):
        # A bitmask to represent the target subset Y.
        # If the j-th bit is 1, a_{j+1} is in Y.
        target_indices = {j for j in range(k) if (i >> j) & 1}
        target_subset = {points[j] for j in target_indices}

        # The construction rule for the formula:
        # The conjunction includes predicates P_j where j is NOT an index
        # of a point in the target subset.
        formula_indices = set(range(k)) - target_indices

        if not formula_indices:
            formula_str = "T (True)"
        else:
            formula_str = " AND ".join([predicates[j] + "(x)" for j in sorted(list(formula_indices))])

        print(f"\nTarget subset Y = {sorted(list(target_subset)) or '{}'}")
        print(f"Proposed formula phi(x) = {formula_str}")

        # Verify the formula for all points
        mismatched = False
        for point_idx in range(k):
            # Evaluate phi(a_{point_idx+1})
            # This is True if P_j(a_{point_idx+1}) is True for all j in formula_indices
            is_phi_true = all(model[point_idx][j] for j in formula_indices)

            # Check if the classification matches the target
            # It should be True if the point is in the target subset, False otherwise
            should_be_true = point_idx in target_indices

            if is_phi_true != should_be_true:
                print(f"  [FAIL] For point {points[point_idx]}: formula evaluates to {is_phi_true}, but should be {should_be_true}")
                mismatched = True
                shattered = False

        if not mismatched:
            print("  [OK] Formula correctly classifies all points.")

    print("\n" + "="*30)
    if shattered:
        print("Conclusion: The set of 4 points was successfully shattered.")
        print("This proves that the VC dimension is at least 4.")
    else:
        print("Conclusion: The set could not be shattered by the construction.")

    print("\nCombining with the proof that VC dim <= k, the final result is:")
    final_k = 4
    print(f"VC Dimension = number of unary predicates = {final_k}")

demonstrate_shattering()