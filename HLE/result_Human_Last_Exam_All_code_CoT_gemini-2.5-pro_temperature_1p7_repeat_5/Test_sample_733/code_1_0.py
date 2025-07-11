import sympy
from sympy.combinatorics.named_groups import (
    CyclicGroup,
    DihedralGroup,
    SymmetricGroup,
    AlternatingGroup,
)
from itertools import combinations

def solve():
    """
    Finds and counts the number of finite groups with maximal by inclusion
    product-free sets of size 2.
    """

    # Helper function to check if a set is product-free
    def is_product_free(s, mult_op):
        """Checks if a subset `s` is product-free."""
        products = {mult_op(x, y) for x in s for y in s}
        return s.isdisjoint(products)

    # Helper function to check if a product-free set is maximal
    def is_maximal_product_free(s, group_elements, mult_op):
        """
        Checks if a product-free set `s` is maximal by inclusion.
        Assumes `s` is already known to be product-free.
        """
        group_elements_set = set(group_elements)
        elements_outside_s = group_elements_set - s
        
        for g in elements_outside_s:
            s_extended = s.union({g})
            if is_product_free(s_extended, mult_op):
                # Found an extension that is still product-free, so s is not maximal.
                return False
                
        # If no extension is product-free, s is maximal.
        return True

    # List of small, non-isomorphic groups to check.
    # This list covers all groups up to order 12 for which sympy has easy constructors,
    # and includes all known groups that satisfy the property.
    groups_to_check = [
        ("C1", CyclicGroup(1)),
        ("C2", CyclicGroup(2)),
        ("C3", CyclicGroup(3)),
        ("C4", CyclicGroup(4)),
        ("C2xC2", DihedralGroup(2)), # Klein four-group
        ("C5", CyclicGroup(5)),
        ("S3", SymmetricGroup(3)), # Also known as D6
        ("C6", CyclicGroup(6)),
        ("C7", CyclicGroup(7)),
        ("C8", CyclicGroup(8)),
        ("D8", DihedralGroup(4)),
        ("A4", AlternatingGroup(4)),
        ("C9", CyclicGroup(9)),
        ("C10", CyclicGroup(10)),
        ("D10", DihedralGroup(5)),
        ("C11", CyclicGroup(11)),
        ("C12", CyclicGroup(12)),
    ]

    found_groups_names = set()

    for name, G in groups_to_check:
        # A group must have at least 3 elements to potentially host a mpfs of size 2,
        # because S must be disjoint from S*S, and e cannot be in S.
        if G.order() < 3:
            continue
            
        elements = list(G.elements)
        mult_op = lambda g1, g2: g1 * g2
        
        # Iterate over all subsets of size 2
        for s_tuple in combinations(elements, 2):
            s = set(s_tuple)
            
            # Check if product-free
            if is_product_free(s, mult_op):
                # If product-free, check if it's maximal
                if is_maximal_product_free(s, elements, mult_op):
                    found_groups_names.add(name)
                    # This group has the property, no need to check other subsets.
                    break 

    print("The script searched for finite groups containing maximal by inclusion product-free sets of size 2.")
    print("The groups found are:")
    # Sort for consistent output
    sorted_names = sorted(list(found_groups_names), key=lambda x: (int(x[1:].split('x')[0]), x))
    for g_name in sorted_names:
        print(f"- {g_name}")

    count = len(found_groups_names)
    print("\nThe final count expressed as a sum is:")
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"{equation_str} = {count}")
    else:
        print("0 = 0")

solve()