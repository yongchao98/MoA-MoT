import itertools

def check_shattering():
    """
    This script verifies that the VC dimension of monotone conjunctions over 4
    variables is at least 4, by showing that a specific set of 4 points
    can be shattered.

    Problem Mapping:
    - n = 4 unary predicates corresponds to n = 4 boolean variables.
    - A point in the universe is characterized by the truth values of the
      4 predicates, so points can be represented as vectors in {0,1}^4.
    - A formula from the logic fragment corresponds to a monotone conjunction
      hypothesis h(x) = x_j1 AND x_j2 AND ...

    The set of points we will shatter, S, consists of 4 vectors in {0,1}^4,
    where each vector has exactly one '0'.
    S = {f_0, f_1, f_2, f_3}
    f_0 = (0, 1, 1, 1)
    f_1 = (1, 0, 1, 1)
    f_2 = (1, 1, 0, 1)
    f_3 = (1, 1, 1, 0)
    """
    n = 4
    
    # Define the set of 4 points to be shattered.
    # We use a list of tuples, where each tuple is a point.
    points = []
    for i in range(n):
        point = [1] * n
        point[i] = 0
        points.append(tuple(point))
        
    print("Set of points S to be shattered:")
    for i, p in enumerate(points):
        print(f"  p{i}: {p}")
    print("-" * 30)

    # We need to show that for every subset of S, there is a hypothesis
    # (a monotone conjunction) that classifies points in the subset as 1
    # and points outside the subset as 0.

    # Iterate through all 2^n = 16 possible subsets of S.
    point_indices = list(range(n))
    total_subsets = 0
    shattered = True
    
    for r in range(n + 1):
        for subset_indices in itertools.combinations(point_indices, r):
            total_subsets += 1
            target_subset_points = {points[i] for i in subset_indices}
            
            # For a target subset of points {p_k | k in K}, the required
            # hypothesis is the conjunction of variables x_j where j is NOT in K.
            # J = K^c
            hypothesis_variable_indices = set(point_indices) - set(subset_indices)

            print(f"Target subset: {{p{k} for k in {set(subset_indices)}}}")
            print(f"  -> Required Hypothesis: Conjunction of variables with indices in {hypothesis_variable_indices}")

            # Verify this hypothesis for all points in S
            mismatch_found = False
            for i, p in enumerate(points):
                # Calculate the hypothesis output for point p
                # h(p) = 1 iff all variables in the conjunction are 1 for p.
                h_output = 1
                for j in hypothesis_variable_indices:
                    if p[j] == 0:
                        h_output = 0
                        break
                
                # Check if the output matches the target classification
                # Target is 1 if p is in the subset, 0 otherwise
                target_output = 1 if p in target_subset_points else 0
                
                if h_output != target_output:
                    print(f"  [FAIL] For p{i}={p}: h(p)={h_output}, but expected {target_output}")
                    mismatch_found = True
                    shattered = False
                    break # Stop checking this hypothesis
            
            if not mismatch_found:
                print("  [SUCCESS] Hypothesis correctly classifies all points in S.")
            
            if not shattered:
                break # Stop checking subsets
        if not shattered:
            break # Stop iterating
            
    print("-" * 30)
    if shattered:
        print(f"All {total_subsets} subsets were successfully generated.")
        print("This proves that the set S is shattered.")
        print("Therefore, the VC dimension is at least 4.")
    else:
        print("Verification failed. The set was not shattered by the implemented logic.")

    print("\nSince it is a known result that the VC dimension of monotone conjunctions")
    print("on n variables is exactly n, and we have shown it is at least 4,")
    print("we conclude that the VC dimension for n=4 is 4.")
    final_equation = "VC_dimension = 4"
    print(f"The final answer is: {final_equation}")


check_shattering()
<<<4>>>