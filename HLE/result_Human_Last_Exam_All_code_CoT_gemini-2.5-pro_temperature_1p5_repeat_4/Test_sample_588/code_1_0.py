import itertools

def find_manifold_tuple():
    """
    This function finds the lexicographically least tuple (a_1, b_1, ..., a_l, b_l)
    with l minimal, such that M(a_i, b_i) is not "full" but their connected sum is.

    A compact, orientable 4-manifold is "full" if its Euler characteristic is 0.
    The problem leads to the Diophantine equation:
      2 * Sum_{i=1 to l} (1-a_i)(1-b_i) = l-1
    where a_i, b_i are non-negative integers not equal to 1.

    This equation requires l to be odd. l=1 is impossible as it forces M(a_1,b_1)
    to be full. The minimal l is 3, giving the equation:
      (1-a_1)(1-b_1) + (1-a_2)(1-b_2) + (1-a_3)(1-b_3) = 1
    
    This script searches for the lexicographically smallest tuple satisfying this.
    """
    l = 3
    target_sum = (l - 1) // 2

    # We search within a small range for the genus values a and b, as we expect
    # the lexicographically smallest solution to use small integers.
    search_range = 5
    candidate_pairs = []
    # Generate pairs (a,b) with a<=b to avoid duplicates, and where a,b != 1.
    for a in range(search_range):
        for b in range(a, search_range):
            if a == 1 or b == 1:
                continue
            # Calculate the corresponding term in the sum
            x_val = (1 - a) * (1 - b)
            candidate_pairs.append(((a, b), x_val))

    # candidate_pairs is now sorted lexicographically, e.g., [((0,0),1), ((0,2),-1), ...]

    # We use itertools.combinations_with_replacement to find three pairs from our
    # candidate list that satisfy the sum. This method yields combinations in
    # lexicographical order, so the first one we find will be the solution.
    for p1, p2, p3 in itertools.combinations_with_replacement(candidate_pairs, 3):
        pair1_val, pair2_val, pair3_val = p1[1], p2[1], p3[1]

        if pair1_val + pair2_val + pair3_val == target_sum:
            # We've found the combination of pairs that works.
            pairs = [p1[0], p2[0], p3[0]]
            
            # The following prints demonstrate that the solution is correct.
            print(f"Found a solution for minimal l = {l}.")
            print(f"The three pairs (a,b) are: {pairs[0]}, {pairs[1]}, and {pairs[2]}.")
            print("Verifying the equation Sum[ (1-a_i)(1-b_i) ] = 1:")
            
            # Print each term of the final equation
            term1 = f"(1-{pairs[0][0]})*(1-{pairs[0][1]})"
            term2 = f"(1-{pairs[1][0]})*(1-{pairs[1][1]})"
            term3 = f"(1-{pairs[2][0]})*(1-{pairs[2][1]})"
            
            print(f"  {term1} = {pair1_val}")
            print(f"  {term2} = {pair2_val}")
            print(f"  {term3} = {pair3_val}")
            print(f"Sum = {pair1_val} + {pair2_val} + {pair3_val} = {target_sum}, which is correct.")

            # Flatten the sorted list of pairs to get the final tuple.
            final_tuple = tuple(item for pair in pairs for item in pair)
            
            print("\nThe lexicographically least tuple is:")
            print(str(final_tuple).replace(" ", ""))
            
            return

find_manifold_tuple()