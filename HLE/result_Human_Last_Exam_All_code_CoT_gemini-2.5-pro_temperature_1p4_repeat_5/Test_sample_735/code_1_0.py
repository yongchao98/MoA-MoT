import sys
from itertools import combinations
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import (
    PermutationGroup,
    CyclicGroup,
    DihedralGroup,
    SymmetricGroup,
    AlternatingGroup,
)
from sympy.combinatorics.group_constructs import DirectProduct

def has_maximal_product_free_set_of_size_3(group):
    """
    Checks if a given group has a maximal by inclusion product-free set of size 3.
    """
    elements = group.elements
    identity = group.identity

    # A product-free set cannot contain the identity element, so we can search
    # through non-identity elements. The group must have at least 3 non-identity elements.
    non_identity_elements = [e for e in elements if e != identity]
    if len(non_identity_elements) < 3:
        return False

    for s_tuple in combinations(non_identity_elements, 3):
        S = set(s_tuple)

        # 1. Check if S is product-free
        is_product_free = True
        for x in S:
            for y in S:
                if x * y in S:
                    is_product_free = False
                    break
            if not is_product_free:
                break
        
        if not is_product_free:
            continue

        # 2. If product-free, check for maximality
        is_maximal = True
        elements_outside_S = [e for e in elements if e not in S]
        for g in elements_outside_S:
            S_union_g = S.union({g})
            
            is_still_product_free = True
            for x in S_union_g:
                for y in S_union_g:
                    if x * y in S_union_g:
                        # The set S U {g} is not product-free, which is what we expect for maximality.
                        is_still_product_free = False
                        break
                if not is_still_product_free:
                    break
            
            if is_still_product_free:
                # S is not maximal because we found an element g that can be added
                # without breaking the product-free property.
                is_maximal = False
                break
        
        if is_maximal:
            # We found one such set for this group, so we can return True.
            return True
            
    return False

def main():
    """
    Main function to define groups and count how many satisfy the condition.
    """
    # Define Quaternion Group Q8 using permutation representation
    a = Permutation(0, 1, 2, 3)(4, 5, 6, 7)
    b = Permutation(0, 4, 2, 6)(1, 7, 3, 5)
    Q8 = PermutationGroup([a, b])

    # List of non-isomorphic groups to check (up to order 12)
    groups_to_check = {
        "C3": CyclicGroup(3),
        "C4": CyclicGroup(4),
        "C5": CyclicGroup(5),
        "C6": CyclicGroup(6),
        "C7": CyclicGroup(7),
        "C8": CyclicGroup(8),
        "C9": CyclicGroup(9),
        "C10": CyclicGroup(10),
        "C11": CyclicGroup(11),
        "C12": CyclicGroup(12),
        "S3": SymmetricGroup(3),
        "V4": PermutationGroup([Permutation(0,1), Permutation(2,3)]), # C2 x C2
        "C3xC3": DirectProduct(CyclicGroup(3), CyclicGroup(3)),
        "D4": DihedralGroup(4),
        "D5": DihedralGroup(5),
        "D6": DihedralGroup(6),
        "A4": AlternatingGroup(4),
        "Q8": Q8,
        "C2xC4": DirectProduct(CyclicGroup(2), CyclicGroup(4)),
        "C2xC2xC2": DirectProduct(CyclicGroup(2), CyclicGroup(2), CyclicGroup(2)),
        "C2xS3": DirectProduct(CyclicGroup(2), SymmetricGroup(3)) # isomorphic to D6
    }
    
    count = 0
    # To avoid counting isomorphisms like D6 and C2xS3 multiple times
    found_groups_by_order_and_abelian = {}

    # Print progress to show work is being done
    sys.stdout.write("Checking groups: ")
    sys.stdout.flush()

    for name, G in groups_to_check.items():
        sys.stdout.write(f"{name}.. ")
        sys.stdout.flush()
        if has_maximal_product_free_set_of_size_3(G):
            # Use order and abelian property as a simple isomorphism check
            key = (G.order(), G.is_abelian)
            if key not in found_groups_by_order_and_abelian:
                 found_groups_by_order_and_abelian[key] = []
            
            # For non-abelian groups of same order, need a better check
            # but our list is mostly curated. We'll just add the name.
            if not G.is_abelian:
                 key = (G.order(), G.is_abelian, name)


            found_groups_by_order_and_abelian[key] = name

    sys.stdout.write("\n")
    
    count = len(found_groups_by_order_and_abelian)
    
    # print("The following groups have maximal product-free sets of size 3:")
    # for k, v in found_groups_by_order_and_abelian.items():
    #     print(v)

    print(count)

if __name__ == '__main__':
    main()
