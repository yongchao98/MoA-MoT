import itertools

def solve_vc_dimension():
    """
    Calculates and demonstrates the VC dimension for the given logic.

    The logic is FO[exist, and, top, bot] over a schema with 4 unary predicates.
    This corresponds to the hypothesis class of positive (monotone) monomials
    over 4 boolean variables, plus the empty concept. The VC dimension is 4.

    This function demonstrates the lower bound (VC >= 4) by showing that a
    carefully chosen set of 4 points can be shattered.
    """

    k = 4
    print(f"Let k be the number of unary predicates. In this problem, k = {k}.")
    print("The hypothesis class corresponds to positive monomials over k variables.")
    print(f"The VC dimension of this class is k.")
    print(f"Therefore, the VC dimension is {k}.\n")
    
    print("--- DEMONSTRATION OF SHATTERING (VC >= 4) ---")
    
    # In this code, we use 0-based indexing for predicates: P_0, P_1, P_2, P_3
    # and points: v_0, v_1, v_2, v_3.
    # The point v_i has a 0 only at position i.
    points = []
    for i in range(k):
        point = [1] * k
        point[i] = 0
        points.append(tuple(point))

    print(f"We will shatter the following set of {k} points:")
    for i, p in enumerate(points):
        print(f"  v{i}: {p}  (Corresponds to properties for predicates P0..P3)")
    print("")

    # Generate all possible 2^k labelings for the k points.
    labels_to_generate = list(itertools.product([1, 0], repeat=k))

    num_shattered = 0
    for i, target_labels in enumerate(labels_to_generate):
        print(f"--- Testing Labeling {i+1}/{2**k} on (v0, v1, v2, v3): {target_labels} ---")
        
        # The all-zeros labeling is achieved by the formula BOT (false).
        if all(label == 0 for label in target_labels):
            print("  This labeling is achieved by the formula: ⊥ (BOT)")
            num_shattered += 1
            continue

        # For any other labeling, we construct the required monomial.
        # Rule: The monomial is a conjunction of predicates P_j for each point v_j
        # that should be labeled 0.
        predicate_indices_in_monomial = []
        for j in range(k):
            if target_labels[j] == 0:
                predicate_indices_in_monomial.append(j)

        # Build the formula string for printing.
        if not predicate_indices_in_monomial:
            # An empty conjunction is always true.
            formula_str = "⊤ (TOP)" 
        else:
            # Note: We use P_j instead of P_j(x) for brevity in the final equation.
            formula_str = " ∧ ".join([f"P_{j}(x)" for j in predicate_indices_in_monomial])
        
        print(f"  Proposed formula for the positive class: {formula_str}")
        
        # Verify that this formula correctly classifies all points.
        actual_labels = []
        for j, point in enumerate(points):
            # A point is classified as 1 (positive) if it satisfies all predicates in the monomial.
            is_positive = 1
            for pred_idx in predicate_indices_in_monomial:
                if point[pred_idx] == 0:
                    is_positive = 0
                    break
            actual_labels.append(is_positive)
        
        print(f"  Applying this formula to points v0..v3 gives labeling: {tuple(actual_labels)}")
        if tuple(actual_labels) == tuple(target_labels):
            print("  Success: The formula produces the target labeling.")
            num_shattered += 1
        else:
            print("  FATAL ERROR: The proof logic is flawed for these points.")
            break
        print("")

    if num_shattered == 2**k:
        print("---------------------------------------------------------------")
        print(f"Successfully generated all {2**k} possible labelings for the set of {k} points.")
        print("This proves the set is shattered, establishing that VC dimension >= 4.")
    
    print("\n--- CONCLUSION ---")
    print(f"Upper Bound (VC <= 4): The hypothesis space has size {2**k + 1}. To shatter d points, we need 2^d <= {2**k + 1}, so d <= 4.")
    print("Lower Bound (VC >= 4): We just demonstrated shattering a set of 4 points.")
    print("Since VC >= 4 and VC <= 4, the VC dimension is exactly 4.")

solve_vc_dimension()