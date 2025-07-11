def multiply_permutations(p1, p2):
    """
    Calculates the composition of two permutations p1 o p2.
    Permutations are represented as tuples showing the 1-indexed mapping.
    e.g., p=(2,1,3) means 1->2, 2->1, 3->3.
    """
    n = len(p1)
    # The new permutation maps i to p1[p2[i-1]-1]
    result = [0] * n
    for i in range(n):
        # p2 is applied first, then p1
        result[i] = p1[p2[i]-1]
    return tuple(result)

def power(p, n):
    """
    Computes the n-th power of a permutation p.
    Handles positive, negative, and zero exponents.
    """
    if n == 0:
        return tuple(range(1, len(p) + 1)) # Identity permutation
    if n < 0:
        # Calculate the inverse permutation
        inv_p = [0] * len(p)
        for i, val in enumerate(p):
            inv_p[val-1] = i+1
        return power(tuple(inv_p), -n)
    
    res = tuple(range(1, len(p) + 1))
    base = p
    while n > 0:
        if n % 2 == 1:
            res = multiply_permutations(base, res)
        base = multiply_permutations(base, base)
        n //= 2
    return res

def get_permutation_cycles(p):
    """
    Finds the cycle decomposition of a permutation p.
    """
    n = len(p)
    not_visited = list(range(1, n + 1))
    cycles = []
    while not_visited:
        cycle = []
        start = not_visited[0]
        current = start
        while current in not_visited:
            cycle.append(current)
            not_visited.remove(current)
            current = p[current - 1]
        cycles.append(tuple(cycle))
    return cycles

def solve_knot_problem():
    """
    Solves the braid problem by calculating the permutation, identifying components,
    and analyzing the relevant sub-braid to determine the knot type.
    """
    n_strands = 5
    braid_word_str = "sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1"
    
    # Define permutations for the generators of the Braid Group B_5
    s1 = (2, 1, 3, 4, 5)
    s2 = (1, 3, 2, 4, 5)
    s3 = (1, 2, 4, 3, 5)
    s4 = (1, 2, 3, 5, 4)

    print(f"Analyzing the braid beta = {braid_word_str} in B_5.")

    # Step 1: Calculate the permutation of the full braid
    # The total permutation is pi = pi(s4^-1) o pi(s3) o pi(s2^2) o pi(s1^2)
    p_s1_sq = power(s1, 2)
    p_s2_sq = power(s2, 2)
    p_s3 = power(s3, 1)
    p_s4_inv = power(s4, -1)
    
    # Multiply permutations in the correct order (top of braid first)
    total_perm = p_s1_sq
    total_perm = multiply_permutations(p_s2_sq, total_perm)
    total_perm = multiply_permutations(p_s3, total_perm)
    total_perm = multiply_permutations(p_s4_inv, total_perm)

    # Step 2: Identify the link components from the permutation cycles
    cycles = get_permutation_cycles(total_perm)
    
    print(f"\nThe permutation of the braid is pi = {cycles}.")
    print(f"The closure of the braid results in a {len(cycles)}-component link.")
    
    unknotted_components = [c for c in cycles if len(c) == 1]
    knotted_component_strands = [c for c in cycles if len(c) > 1][0]
    
    print(f"The link components correspond to the cycles: {unknotted_components[0]}, {unknotted_components[1]}, and {knotted_component_strands}.")
    print("As given, two components are unknots, corresponding to the single-strand cycles.")

    # Step 3: Analyze the sub-braid for the knotted component
    print("\nThe component of interest is the one formed by strands {3, 4, 5}.")
    print("Its intrinsic knot type is determined by the sub-braid that acts only on these strands.")
    print("This sub-braid is sigma_3 * sigma_4^-1.")
    
    # Step 4: Identify the knot type
    print("\nTo identify the knot, we can analyze the equivalent 3-strand braid: sigma_1 * sigma_2^-1.")
    sub_braid_exponent_sum = 1 + (-1)
    print(f"The exponent sum of this sub-braid is 1 + (-1) = {sub_braid_exponent_sum}.")
    print("An exponent sum of 0 rules out chiral knots like the Trefoil, whose exponent sum must be odd and non-zero. This suggests the unknot or an amphichiral knot like the Figure-8.")
    
    print("\nThe braid sigma_1 * sigma_2^-1 is a well-known representation of the Figure-8 knot (4_1).")
    print("\nTherefore, the third connected component is equivalent to the Figure-8 knot.")

# Run the analysis
solve_knot_problem()