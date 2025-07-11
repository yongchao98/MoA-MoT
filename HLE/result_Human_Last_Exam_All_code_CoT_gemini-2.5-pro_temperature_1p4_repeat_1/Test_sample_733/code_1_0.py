import sympy.combinatorics.perm_groups as perm
from itertools import combinations

def has_maximal_product_free_set_of_size_2(G):
    """
    Checks if a group G has a maximal by inclusion product-free set of size 2.
    """
    # Trivial groups or groups of order 2 cannot have a product-free set of size 2.
    if G.order() < 3:
        return False
    
    elements = list(G.elements)
    identity = G.identity

    # Iterate over all subsets of size 2
    for s_tuple in combinations(elements, 2):
        S = set(s_tuple)
        
        # --- Check if S is product-free ---
        is_pf = True
        # A product-free set cannot contain the identity element, because if e is in S, e*e=e is in S.
        if identity in S:
            continue
            
        # Check all products S·S
        for x in S:
            for y in S:
                if (x * y) in S:
                    is_pf = False
                    break
            if not is_pf:
                break
        
        if not is_pf:
            continue
        
        # --- If product-free, check for maximality ---
        is_maximal = True
        # For any g not in S, S U {g} must NOT be product-free.
        elements_outside_S = [g for g in elements if g not in S]
        for g in elements_outside_S:
            S_prime = S.union({g})
            
            # Check if S_prime is product-free
            is_s_prime_pf = True
            for x in S_prime:
                for y in S_prime:
                    if (x * y) in S_prime:
                        is_s_prime_pf = False
                        break
                if not is_s_prime_pf:
                    break

            # If S_prime IS product-free, then S is NOT maximal.
            if is_s_prime_pf:
                is_maximal = False
                break
        
        if is_maximal:
            return True # Found one such set, so the group qualifies.
            
    return False

def find_groups():
    """
    Tests a library of small groups and prints the ones that have a
    maximal product-free set of size 2.
    """
    # A dictionary of group names and their SymPy permutation group representations.
    # This covers all groups up to order 12.
    groups_to_check = {
        "C2": perm.CyclicGroup(2),
        "C3": perm.CyclicGroup(3),
        "C4": perm.CyclicGroup(4),
        "C2xC2": perm.AbelianGroup(2,2),
        "C5": perm.CyclicGroup(5),
        "C6": perm.CyclicGroup(6),
        "S3 (or D6)": perm.SymmetricGroup(3),
        "C7": perm.CyclicGroup(7),
        "C8": perm.CyclicGroup(8),
        "C4xC2": perm.AbelianGroup(4,2),
        "C2xC2xC2": perm.AbelianGroup(2,2,2),
        "D8": perm.DihedralGroup(4),
        "Q8": perm.QuaternionGroup(),
        "C9": perm.CyclicGroup(9),
        "C3xC3": perm.AbelianGroup(3,3),
        "C10": perm.CyclicGroup(10),
        "D10": perm.DihedralGroup(5),
        "C11": perm.CyclicGroup(11),
        "C12": perm.CyclicGroup(12),
        "C6xC2": perm.AbelianGroup(6,2),
        "A4": perm.AlternatingGroup(4),
        "D12": perm.DihedralGroup(6),
        "Q12": perm.DicyclicGroup(3)
    }

    found_groups_count = 0
    print("Searching for finite groups with maximal product-free sets of size 2...")
    
    for name, G in groups_to_check.items():
        if has_maximal_product_free_set_of_size_2(G):
            print(f"Found one: {name} (order {G.order()})")
            found_groups_count += 1
    
    print("\nTotal number of such groups found:")
    print(found_groups_count)

if __name__ == '__main__':
    find_groups()
