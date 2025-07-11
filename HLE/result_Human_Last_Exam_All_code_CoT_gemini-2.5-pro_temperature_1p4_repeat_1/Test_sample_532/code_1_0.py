import itertools

def find_mips_and_union_for_cyclic_group(n):
    """
    Finds all Maximal by Inclusion Product-Free Sets (MIPS) for a cyclic
    group C_n and computes their union.

    Args:
        n (int): The order of the cyclic group C_n = {0, 1, ..., n-1}.
    """
    if n <= 1:
        print(f"The group C_{n} is trivial or ill-defined for this problem.")
        return

    elements = set(range(n))
    print(f"Analyzing the finite nilpotent group C_{n} = {elements}")
    print(f"The identity element is 0.")
    print("-" * 30)

    # --- Step 1: Find all product-free sets ---
    # We generate all non-empty subsets of the group elements.
    all_subsets = []
    # Note: frozenset is used so sets can be added to other sets.
    for r in range(1, n + 1):
        for s_tuple in itertools.combinations(elements, r):
            all_subsets.append(frozenset(s_tuple))

    product_free_sets = []
    for s in all_subsets:
        is_pf = True
        # A set S is product-free if for all x, y in S, x+y is not in S.
        for x in s:
            for y in s:
                if (x + y) % n in s:
                    is_pf = False
                    break
            if not is_pf:
                break
        if is_pf:
            product_free_sets.append(s)

    # --- Step 2: Filter for Maximal product-free sets (MIPS) ---
    mips = []
    for s1 in product_free_sets:
        is_maximal = True
        for s2 in product_free_sets:
            # Check if s1 is a proper subset of any other product-free set s2
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            mips.append(s1)

    print("The Maximal by Inclusion Product-Free Sets (MIPS) are:")
    for s in sorted(list(mips), key=len):
        print(f"  {sorted(list(s))}")

    # --- Step 3: Compute and print the union of the MIPS ---
    union_of_mips = set().union(*mips)

    print("-" * 30)
    print("The final union of all the MIPS is:")
    # The prompt asked to "output each number in the final equation"
    # which we interpret as printing the elements of the resulting set.
    print(f"U(MIPS) = {sorted(list(union_of_mips))}")

    # --- Step 4: Final verification ---
    non_identity_elements = elements - {0}
    if union_of_mips == non_identity_elements:
        print("\nConclusion: The union covers all elements except the identity, as expected.")
    else:
        print("\nConclusion: The union does NOT cover all non-identity elements.")

# Demonstrate with the cyclic group C_6, which is a finite nilpotent group.
find_mips_and_union_for_cyclic_group(6)
