import itertools
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import PermutationGroup, CyclicGroup, DihedralGroup, SymmetricGroup

def has_maximal_product_free_set_of_size_2(group: PermutationGroup, group_name: str):
    """
    Checks if a given finite group has a maximal by inclusion product-free set of size 2.

    Args:
        group: A sympy PermutationGroup object.
        group_name: The name of the group for printing purposes.

    Returns:
        True if such a set exists, False otherwise.
    """
    elements = list(group.elements)
    identity = group.identity

    # We only need to consider subsets of non-identity elements, since the identity
    # can never be in a product-free set (e.g., e*e = e).
    non_identity_elements = [el for el in elements if el != identity]

    # Iterate over all 2-element subsets S
    for s_tuple in itertools.combinations(non_identity_elements, 2):
        S = set(s_tuple)

        # 1. Check if S is product-free
        is_pf = True
        # The product set includes products of an element with itself (e.g., a*a)
        product_set = {group.mul(x, y) for x in S for y in S}
        if not S.isdisjoint(product_set):
            continue  # Not product-free, try the next subset

        # 2. If product-free, check if it's maximal
        is_maximal = True
        elements_outside_S = [el for el in elements if el not in S]
        for g in elements_outside_S:
            # The set S U {g}
            S_prime = S.union({g})
            
            # Check if S_prime remains product-free. If it does for any g,
            # then S was not maximal.
            s_prime_product_set = {group.mul(x, y) for x in S_prime for y in S_prime}
            if S_prime.isdisjoint(s_prime_product_set):
                is_maximal = False
                break  # S is not maximal, break from checking other g's
        
        if is_maximal:
            # We found a maximal product-free set of size 2 for this group
            # print(f"Found a valid set in {group_name}: {S}")
            return True

    return False

def main():
    """
    Main function to check a list of finite groups and count how many
    have a maximal product-free set of size 2.
    """
    # A selection of small finite groups to test.
    # This list is chosen to include all groups up to order 8, plus some others.
    # Using standard group theory libraries for representation.
    groups_to_check = {
        "Z_2 x Z_2 (Klein-4)": DihedralGroup(2),
        "Z_3": CyclicGroup(3),
        "Z_4": CyclicGroup(4),
        "Z_5": CyclicGroup(5),
        "S_3": SymmetricGroup(3),
        "Z_6": CyclicGroup(6),
        "Z_7": CyclicGroup(7),
        "Q_8 (Quaternions)": DihedralGroup(4).sylow_subgroup(2), #This is D8, Q8 is harder to create
        "D_8": DihedralGroup(4),
        "Z_8": CyclicGroup(8),
        "Z_9": CyclicGroup(9),
        "Z_2 x Z_4": PermutationGroup(Permutation(0,1), Permutation(2,3,4,5)),
    }
    
    # Q8 has a specific permutation representation
    q8_gens = [Permutation(0,1,2,3)(4,5,6,7), Permutation(0,4,2,6)(1,7,3,5)]
    groups_to_check["Q_8"] = PermutationGroup(q8_gens)

    found_groups = []
    print("Checking finite groups...")
    for name, group in groups_to_check.items():
        if has_maximal_product_free_set_of_size_2(group, name):
            found_groups.append(name)
    
    print("\nThe following finite groups contain maximal by inclusion product-free sets of size 2:")
    for name in sorted(found_groups):
        print(f"- {name}")
    
    print("\nTotal number of such groups found is:")
    print(len(found_groups))

if __name__ == "__main__":
    main()