import itertools

def demonstrate_s3_not_filled():
    """
    This function demonstrates that the group S_3 is not a filled group.
    S_3 is the nonabelian group of order 6, which fits the form 2*q^m with q=3, m=1.
    """
    # Step 1: Define the group S_3
    # Elements are represented by indices 0-5 for computation.
    # e, (12), (13), (23), (123), (132)
    elements = ['e', '(12)', '(13)', '(23)', '(123)', '(132)']
    elem_to_idx = {name: i for i, name in enumerate(elements)}
    G_indices = set(range(6))

    # Cayley table (multiplication table) for S_3
    # cayley[i][j] = i * j
    cayley = [
        [0, 1, 2, 3, 4, 5],
        [1, 0, 5, 4, 3, 2],
        [2, 4, 0, 5, 1, 3],
        [3, 5, 4, 0, 2, 1],
        [4, 2, 3, 1, 5, 0],
        [5, 3, 1, 2, 0, 4]
    ]
    
    # Inverse of each element
    inverses = [0, 1, 2, 3, 5, 4]

    print("--- Analyzing the group S_3 (order 6) ---")
    print(f"Group G = { {elements[i] for i in G_indices} }")
    print(f"Identity element e = '{elements[0]}'")

    # Step 2: Find a candidate maximal product-free set S.
    # A known candidate is the set of all elements of order 2 (involutions).
    S_names = {'(12)', '(13)', '(23)'}
    S = {elem_to_idx[name] for name in S_names}
    print(f"\nLet's test the set S = {S_names}, which contains all elements of order 2.")

    # Step 3: Verify that S is product-free.
    # A set S is product-free if for all x, y in S, x*y is not in S.
    is_product_free = True
    for x in S:
        for y in S:
            product = cayley[x][y]
            if product in S:
                is_product_free = False
                print(f"Verification failed: S is not product-free because {elements[x]} * {elements[y]} = {elements[product]} is in S.")
                break
        if not is_product_free:
            break
    
    if is_product_free:
        print("Verification passed: S is a product-free set.")

    # Step 4: Verify that S is maximal.
    # S is maximal if for any element g not in S, the set S U {g} is not product-free.
    is_maximal = True
    elements_outside_S = G_indices - S - {elem_to_idx['e']} # Exclude identity
    for g in elements_outside_S:
        S_plus_g = S.union({g})
        is_still_product_free = True
        # Check all pairs in the new set
        for x in S_plus_g:
            for y in S_plus_g:
                product = cayley[x][y]
                if product in S_plus_g:
                    # Found a product within the new set, so it's not product-free
                    print(f"Testing maximality with element '{elements[g]}':")
                    print(f"  The set S U {{'{elements[g]}'}} is not product-free because {elements[x]} * {elements[y]} = {elements[product]}.")
                    is_still_product_free = False
                    break
            if not is_still_product_free:
                break
        if is_still_product_free:
            is_maximal = False
            print(f"Verification failed: S is not maximal, as S U {{'{elements[g]}'}} is still product-free.")
            break

    if is_maximal:
        print("Verification passed: S is a maximal product-free set.")

    # Step 5: Check if this maximal product-free set S fills the group G.
    # The "fills" condition is G = S U S_inv U {e}.
    print("\n--- Checking the 'fills' condition: G = S U S⁻¹ U {e} ---")
    
    S_inv = {inverses[i] for i in S}
    S_union_S_inv_union_e = S.union(S_inv).union({elem_to_idx['e']})
    
    G_set_names = {elements[i] for i in G_indices}
    S_set_names = {elements[i] for i in S}
    S_inv_set_names = {elements[i] for i in S_inv}
    e_set_name = {elements[elem_to_idx['e']]}
    Union_set_names = {elements[i] for i in S_union_S_inv_union_e}

    print(f"The group G is: {G_set_names}")
    print(f"The maximal product-free set S is: {S_set_names}")
    print(f"The set of inverses S⁻¹ is: {S_inv_set_names}")
    print(f"The identity set {{e}} is: {e_set_name}")
    print("\nThe final equation we are testing is G = S U S⁻¹ U {e}")
    print(f"The union S U S⁻¹ U {{e}} results in the set: {Union_set_names}")
    
    print(f"\nComparing the sizes:")
    print(f"|G| = {len(G_indices)}")
    print(f"|S U S⁻¹ U {{e}}| = {len(S_union_S_inv_union_e)}")

    if G_indices == S_union_S_inv_union_e:
        print("\nConclusion: The set S fills G.")
    else:
        print("\nConclusion: The set S does NOT fill G because the sets are not equal.")
        print("Since we found a maximal product-free set that does not fill the group, S_3 is not a filled group.")

if __name__ == '__main__':
    demonstrate_s3_not_filled()
