import itertools

def solve_vc_dimension():
    """
    Calculates and explains the VC dimension for the given first-order logic fragment.
    """
    k = 4  # Number of unary predicates in schema S
    predicates = list(range(k))
    points = list(range(k))

    print("Step 1: Analysis of the Logic")
    print("The logic is FO_exists,and,T,F[S] with S having 4 unary predicates (P0, P1, P2, P3).")
    print("Any formula phi(x) in this logic simplifies to either F (False) or a monomial of the form: P_i(x) AND P_j(x) ...")
    print("This is because any existential quantifiers on other variables (e.g., exists y. P_k(y)) result in a sentence that is either True or False in a model,")
    print("either making the whole formula equivalent to the part with x, or making it False.\n")
    print(f"Thus, the problem is to find the VC dimension for concepts defined by monomials over {k} predicates.\n")

    print(f"------------------------------------")
    print(f"Step 2: Show VC dimension >= {k}")
    print(f"------------------------------------")
    print(f"We will show that a set X of {k} points, let's call them c0, c1, c2, c3, can be shattered.")
    print("To do this, we construct a model M where P_i(c_j) is true if and only if i != j.")
    print("Let's visualize this model as a matrix A, where A[i, j] is 1 if P_i(c_j) is true, and 0 otherwise:")
    
    # Define the model M(i, j) which evaluates P_i(c_j)
    def model(predicate_idx, point_idx):
        return predicate_idx != point_idx

    for i in range(k):
        row = [int(model(i, j)) for j in range(k)]
        print(f"  P{i}: {row}")

    print("\nNow we show that for any target subset of X, there is a formula that defines it.")
    
    # Generate all 2^k subsets of points
    subsets_of_points = []
    for i in range(len(points) + 1):
        for subset in itertools.combinations(points, i):
            subsets_of_points.append(set(subset))

    # For each target subset, we find the formula (a conjunction of predicates) that generates it
    print("\n{:<25} {:<25}".format("Target Subset of Points", "Generating Formula"))
    print("-" * 50)

    for target_subset_indices in subsets_of_points:
        # The required set of predicates is the complement of the point indices
        predicate_indices_for_formula = set(predicates) - target_subset_indices
        
        # Build the formula string for display
        if not predicate_indices_for_formula:
            formula_str = "T (True)"
        else:
            formula_str = " AND ".join(sorted([f"P{p}" for p in predicate_indices_for_formula]))
        
        target_subset_str = "{" + ", ".join(sorted([f"c{p}" for p in target_subset_indices])) + "}" if target_subset_indices else "{}"
        print("{:<25} {:<25}".format(target_subset_str, formula_str))
        
    print(f"\nSince all {2**k} subsets of X can be generated, the set is shattered.")
    print(f"This proves that the VC dimension is at least {k}.\n")


    print(f"------------------------------------")
    print(f"Step 3: Show VC dimension < {k + 1}")
    print(f"------------------------------------")
    num_formulas = 2**k
    num_concepts = num_formulas + 1 # +1 for the False formula
    
    d = k + 1
    num_subsets_to_shatter = 2**d
    
    print(f"To shatter a set of d = {d} points, we would need to generate 2^{d} = {num_subsets_to_shatter} distinct subsets.")
    print(f"Our language consists of conjunctions of a subset of {k} predicates.")
    print(f"The number of such subsets of predicates is 2^{k} = {num_formulas}.")
    print(f"This gives {num_formulas} monomial formulas. Including the special formula F (False), our language can define at most {num_concepts} distinct concepts in any model.")
    print(f"To shatter {d} points, we need {num_subsets_to_shatter} concepts, but we only have {num_concepts}.")
    print(f"Since {num_subsets_to_shatter} > {num_concepts}, it's impossible to shatter a set of {d} points.")
    print(f"This proves that the VC dimension is less than {d}.\n")
    
    print(f"------------------------------------")
    print("Step 4: Conclusion")
    print(f"------------------------------------")
    print(f"From Step 2, we have VC-dim >= {k}.")
    print(f"From Step 3, we have VC-dim < {k + 1}.")
    print(f"Therefore, the VC dimension of FO_exists,and,T,F[S] with {k} unary predicates is exactly {k}.")

solve_vc_dimension()