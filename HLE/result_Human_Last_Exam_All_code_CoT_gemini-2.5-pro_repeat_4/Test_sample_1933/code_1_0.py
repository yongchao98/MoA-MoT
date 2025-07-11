import itertools

def solve_vc_dimension():
    """
    This script explains and calculates the VC dimension for the given problem.
    The task is to find the VC dimension of formulas in FO[exists, AND] over a
    schema with 4 unary predicates.
    """

    # The number of unary predicates in the schema S
    num_predicates = 4

    print(f"The schema has {num_predicates} unary predicates (P1, P2, P3, P4).")
    print("The problem is equivalent to finding the VC dimension of concepts defined by monomials (conjunctions of literals) over these 4 predicates.")
    print(f"The VC dimension of this class is known to be equal to the number of predicates, which is {num_predicates}.\n")

    print(f"--- Proof of VC Dimension >= {num_predicates} ---")
    print(f"We prove this by showing that a specific set of {num_predicates} points can be 'shattered'.")
    print("To shatter a set means we can use our formulas to select any possible subset of it.")
    print("\nLet's define a set of 4 points: {u1, u2, u3, u4}.")
    print("Let's define a model where P_i(u_j) is TRUE if and only if i = j.")
    print("This makes the points independent, like basis vectors in a 4D space:")
    print("u1 satisfies (P1=T, P2=F, P3=F, P4=F)")
    print("u2 satisfies (P1=F, P2=T, P3=F, P4=F)")
    print("u3 satisfies (P1=F, P2=F, P3=T, P4=F)")
    print("u4 satisfies (P1=F, P2=F, P4=T, P4=T)\n")


    print(f"We will now show how to generate all 2^{num_predicates} = {2**num_predicates} subsets of {{u1, u2, u3, u4}} using a monomial formula.\n")

    points = [f"u{i+1}" for i in range(num_predicates)]

    # Iterate through all 2^d possible subsets
    num_subsets = 2**num_predicates
    for i in range(num_subsets):
        target_subset_points = []
        target_indices = set()  # Use 1-based indexing for predicates
        for j in range(num_predicates):
            if (i >> j) & 1:
                target_subset_points.append(points[j])
                target_indices.add(j + 1)

        # Format the target subset for printing
        if not target_subset_points:
            target_subset_str = "{}"
        else:
            target_subset_str = "{" + ", ".join(target_subset_points) + "}"

        # The general formula to select a target subset S is: AND_{k where u_k is NOT in S} (NOT P_k(x))
        all_indices = set(range(1, num_predicates + 1))
        indices_for_negated_literals = all_indices - target_indices

        literals = [f"NOT P{k}(x)" for k in sorted(list(indices_for_negated_literals))]

        if not literals:
            formula = "TRUE" # The empty conjunction is TRUE
        else:
            formula = " AND ".join(literals)

        print(f"To select the subset {target_subset_str:<22} use the formula: {formula}")

    print("\nSince every possible subset can be isolated, the set {u1, u2, u3, u4} is shattered.")
    print(f"This proves that the VC dimension is at least {num_predicates}.\n")

    print(f"--- Proof of VC Dimension <= {num_predicates} ---")
    print(f"It's a known result that for monomials over d attributes, no set of d+1 points can be shattered.")
    print(f"In our case, with d={num_predicates}, no set of {num_predicates + 1} points can be shattered.")
    print("This provides the upper bound on the VC dimension.\n")

    print("--- Conclusion ---")
    print(f"Since the VC dimension is >= {num_predicates} and <= {num_predicates}, it must be exactly {num_predicates}.")
    print(f"Final Answer: The VC dimension is {num_predicates}.")

solve_vc_dimension()