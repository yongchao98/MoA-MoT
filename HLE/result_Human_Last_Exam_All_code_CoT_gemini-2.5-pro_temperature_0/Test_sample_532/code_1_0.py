import itertools

def is_product_free(s, table):
    """Checks if a subset s is product-free given a Cayley table."""
    if not s:
        return True
    for e1 in s:
        for e2 in s:
            product = table[e1][e2]
            if product in s:
                return False
    return True

def find_mips_and_check_filled(group_name, order):
    """
    Finds all Maximal by Inclusion Product-Free Sets (MIPS) for the cyclic group C_n
    and checks if the group is filled.
    """
    # Generate the Cayley table for the cyclic group C_n (additive mod n)
    table = [[(i + j) % order for j in range(order)] for i in range(order)]
    
    elements = list(range(order))
    # The identity element is 0. Non-identity elements are 1, 2, ..., order-1.
    non_identity_elements = set(range(1, order))

    print(f"--- Analyzing Group: {group_name} ---")
    print(f"This is the cyclic group of order {order}. Elements are {elements}, with 0 as the identity.")
    print(f"The set of non-identity elements G \\ {{e}} is: {non_identity_elements}\n")

    # 1. Find all non-empty product-free subsets of G \ {e}
    all_pf_sets = []
    # Iterate through all possible non-empty subsets of non-identity elements
    for i in range(1, order):
        for subset_tuple in itertools.combinations(non_identity_elements, i):
            subset = set(subset_tuple)
            if is_product_free(subset, table):
                all_pf_sets.append(subset)

    # 2. From all product-free sets, find the maximal ones (MIPS)
    mips = []
    for s1 in all_pf_sets:
        is_maximal = True
        for s2 in all_pf_sets:
            # Check if s1 is a proper subset of s2
            if s1.issubset(s2) and s1 != s2:
                is_maximal = False
                break
        if is_maximal:
            mips.append(s1)

    print("Found the following Maximal Product-Free Sets (MIPS):")
    if not mips:
        print("None")
    else:
        for s in sorted([sorted(list(m)) for m in mips]):
            print(s)

    # 3. Take the union of all MIPS and check if the group is filled
    union_of_mips = set().union(*mips)
    uncovered_elements = non_identity_elements - union_of_mips

    print(f"\nThe union of all MIPS is: {sorted(list(union_of_mips)) if union_of_mips else '{}'}")
    
    print("\n--- Conclusion ---")
    if not uncovered_elements:
        print(f"The union of MIPS covers all non-identity elements.")
        print(f"Therefore, {group_name} is a filled group.")
    else:
        print(f"The union of MIPS does NOT cover all non-identity elements.")
        print(f"The uncovered elements are: {sorted(list(uncovered_elements))}")
        print(f"Therefore, {group_name} is NOT a filled group.")
    print("-" * 40)


# The mathematical theorem states that a finite nilpotent group is filled
# if and only if it is a 2-group (i.e., its order is a power of 2).

# We demonstrate this with two examples of abelian (and therefore nilpotent) groups.

# Example 1: C_4 (Cyclic group of order 4). Order is 2^2, so it's a 2-group.
# The theorem predicts it IS filled.
find_mips_and_check_filled("C_4", 4)

# Example 2: C_6 (Cyclic group of order 6). Order is not a power of 2.
# The theorem predicts it IS NOT filled.
find_mips_and_check_filled("C_6", 6)
