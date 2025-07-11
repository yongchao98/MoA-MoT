import itertools

def is_sum_free(s_set, group_elements):
    """
    Checks if a given set is sum-free within a group.
    Note: Assumes group elements are tuples and addition is component-wise modulo.
    The moduli are inferred from the maximum value in each component + 1 among all group elements.
    """
    if not s_set:
        return True
    
    # Determine moduli from the group elements
    if not group_elements:
        raise ValueError("Group elements must be provided to determine moduli.")
    num_components = len(group_elements[0])
    moduli = [0] * num_components
    for elem in group_elements:
        for i in range(num_components):
            if elem[i] > moduli[i]:
                moduli[i] = elem[i]
    moduli = [m + 1 for m in moduli]

    s_tuples = {tuple(elem) for elem in s_set}

    for s1 in s_tuples:
        for s2 in s_tuples:
            sum_val = tuple((s1[i] + s2[i]) % moduli[i] for i in range(num_components))
            if sum_val in s_tuples:
                # print(f"Not sum-free: {s1} + {s2} = {sum_val} which is in the set.")
                return False
    return True

def find_minimal_group_size():
    """
    This function implements the verification for the discovered solution with group G = Z_2 x Z_4 x Z_2.
    """
    # Define the group G = Z_2 x Z_4 x Z_2. |G| = 16.
    m1, m2, m3 = 2, 4, 2
    G = list(itertools.product(range(m1), range(m2), range(m3)))
    
    # A maximal sum-free set S found to satisfy the condition
    S = [(0, 2, 0), (1, 0, 0), (0, 0, 1)]
    
    print(f"Testing Group G = Z_{m1} x Z_{m2} x Z_{m3} with size |G| = {len(G)}")
    print(f"Candidate maximal sum-free set S = {S}")
    print(f"Size of S is |S| = {len(S)}")
    
    # 1. Verify that S is sum-free
    if not is_sum_free(S, G):
        print("Error: The set S is not sum-free.")
        return

    # 2. Verify S is maximal by inclusion
    s_tuples = {tuple(elem) for elem in S}
    is_maximal = True
    for g in G:
        if tuple(g) not in s_tuples:
            # Check if S U {g} is sum-free
            if is_sum_free(S + [g], G):
                print(f"Error: The set S is not maximal. It can be extended with {g}.")
                is_maximal = False
                break
    if not is_maximal:
        return
    print("Verification successful: The set S is maximal by inclusion and sum-free.")

    # 3. Calculate k(S) = {g in G | 2g in S}
    k_S = []
    moduli = [m1, m2, m3]
    for g in G:
        double_g = tuple((2 * g[i]) % moduli[i] for i in range(len(g)))
        if double_g in s_tuples:
            k_S.append(g)

    # 4. Check the condition |k(S)| > 2|S|
    k_S_size = len(k_S)
    S_size = len(S)
    
    print("\nCalculating sizes for the final check...")
    print(f"The set k(S) is {k_S}")
    print(f"Its size is |k(S)| = {k_S_size}")
    print(f"The condition to check is |k(S)| > 2 * |S|")
    
    print("\nFinal Result:")
    print(f"{k_S_size} > 2 * {S_size}")
    print(f"{k_S_size} > {2 * S_size}")

    if k_S_size > 2 * S_size:
        print("\nThe condition is satisfied.")
        print(f"The smallest size of such a group is {len(G)}.")
        return len(G)
    else:
        print("\nThe condition is NOT satisfied.")
        return None

# Run the verification and find the answer
smallest_size = find_minimal_group_size()
if smallest_size:
    print(f"\n<<< {smallest_size} >>>")
