import itertools

def solve_manifold_problem():
    """
    This function finds the lexicographically least tuple (a_1, b_1, ..., a_ell, b_ell)
    that satisfies the conditions of the problem by performing a systematic search.
    The logic is explained in the plan above. The minimal number of manifolds is l=3,
    leading to the equation:
    (1-a_1)(1-b_1) + (1-a_2)(1-b_2) + (1-a_3)(1-b_3) = 1
    where a_i, b_i are non-negative integers not equal to 1.
    """

    # We search for the solution with the minimal ell=3.
    # We need to find the lexicographically smallest tuple, so we search for pairs (a,b)
    # in lexicographical order.
    # The search range for a and b can be kept small, as the smallest solution
    # is expected to use small numbers.
    search_values = [0, 2, 3, 4]  # Possible values for a_i, b_i (must not be 1)
    candidate_pairs = []
    for a in search_values:
        for b in search_values:
            pair = (a, b)
            x = (1 - a) * (1 - b)
            candidate_pairs.append({'pair': pair, 'x': x})

    # Use combinations_with_replacement to find three pairs whose x values sum to 1.
    # This iterator yields combinations in sorted order based on the input list.
    # Since `candidate_pairs` is sorted lexicographically by `pair`, the first
    # solution we find will correspond to the lexicographically smallest tuple.
    for p1, p2, p3 in itertools.combinations_with_replacement(candidate_pairs, 3):
        if p1['x'] + p2['x'] + p3['x'] == 1:
            # This is the first, and therefore lexicographically smallest, solution.
            
            # Define the manifolds for the explanation.
            pair1, pair2, pair3 = p1['pair'], p2['pair'], p3['pair']
            M1, M2, M3 = f"M{pair1}", f"M{pair2}", f"M{pair3}"
            
            # Calculate their Euler characteristics.
            chi1 = 4 * p1['x']
            chi2 = 4 * p2['x']
            chi3 = 4 * p3['x']
            
            # The connect-sum is M = M1 # M2 # M3. Its Euler characteristic is:
            # chi(M) = chi(M1) + chi(M2) + chi(M3) - 2*(l-1)
            # For l=3, this is chi(M) = chi1 + chi2 + chi3 - 4.
            # We verify this is 0.
            
            print("The minimal number of manifolds is l=3.")
            print(f"The three manifolds M(a,b) that are not full are {M1}, {M2}, and {M3}.")
            print(f"Their Euler characteristics are chi({M1}) = {chi1}, chi({M2}) = {chi2}, and chi({M3}) = {chi3}.")
            print("The Euler characteristic of their connect-sum M is given by the equation:")
            print(f"chi(M) = chi({M1}) + chi({M2}) + chi({M3}) - 2*(3-1)")
            print(f"chi(M) = {chi1} + {chi2} + {chi3} - 4")
            print(f"chi(M) = {chi1 + chi2 + chi3} - 4 = 0")
            print("The signature of the connect-sum is also 0, so the resulting manifold is full.")
            
            # Construct the final tuple.
            result_tuple = pair1 + pair2 + pair3
            
            print("\nThe lexicographically least tuple is:")
            # Print the tuple in the requested format (flat, no spaces).
            print(str(result_tuple).replace(" ", ""))
            return

solve_manifold_problem()