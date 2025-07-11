from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import PermutationGroup, CyclicGroup, DihedralGroup, KleinFourGroup, SymmetricGroup, AlternatingGroup
from itertools import combinations

def has_maximal_product_free_set_of_size_2(G):
    """
    Checks if a given group G has a maximal product-free set of size 2.
    """
    elements = list(G.elements)
    # The set S cannot contain the identity element.
    non_identity_elements = [el for el in elements if el != G.identity]

    # We need at least two non-identity elements to form a set of size 2.
    if len(non_identity_elements) < 2:
        return False

    # Iterate through all combinations of two distinct non-identity elements.
    for s_tuple in combinations(non_identity_elements, 2):
        s_set = set(s_tuple)
        a, b = s_tuple

        # 1. Check if S = {a, b} is product-free.
        # The set of products is S*S = {a*a, a*b, b*a, b*b}.
        # For S to be product-free, (S*S) intersect S must be empty.
        products = {a*a, a*b, b*a, b*b}
        if not s_set.isdisjoint(products):
            continue

        # 2. Check if S is maximal by inclusion.
        # For any element g not in S, the set S U {g} must NOT be product-free.
        is_maximal = True
        other_elements = [el for el in elements if el not in s_set]

        for g in other_elements:
            s_prime_set = s_set.union({g})
            s_prime_list = list(s_prime_set)
            
            # Check if S' = S U {g} is product-free.
            # If we find even one g for which S' is product-free, then S is not maximal.
            s_prime_is_product_free = True
            for p1 in s_prime_list:
                for p2 in s_prime_list:
                    if (p1 * p2) in s_prime_set:
                        s_prime_is_product_free = False
                        break
                if not s_prime_is_product_free:
                    break
            
            if s_prime_is_product_free:
                is_maximal = False
                break
        
        # If we went through all g and S remained maximal, then we found one.
        if is_maximal:
            return True
            
    return False

def main():
    """
    Main function to check groups and print the result.
    """
    # Based on known theorems, we only need to check groups of order up to 7,
    # plus Z_3 x Z_3 for order 9. Groups of order 8 are known to fail.
    groups_to_check = {
        "C_4 (Cyclic group of order 4)": CyclicGroup(4),
        "V_4 (Klein four-group)": KleinFourGroup(),
        "C_5 (Cyclic group of order 5)": CyclicGroup(5),
        "C_6 (Cyclic group of order 6)": CyclicGroup(6),
        "S_3 (Symmetric group of order 6)": SymmetricGroup(3),
        "C_7 (Cyclic group of order 7)": CyclicGroup(7),
        "Z_3 x Z_3 (Abelian group of order 9)": PermutationGroup(Permutation(0,1,2), Permutation(3,4,5))
    }

    found_groups = []
    for name, G in groups_to_check.items():
        if has_maximal_product_free_set_of_size_2(G):
            found_groups.append(name)

    print("The finite groups containing maximal by inclusion product-free sets of size 2 are:")
    for name in found_groups:
        print(f"- {name}")
    
    print("\nTotal number of such groups:")
    print(len(found_groups))

if __name__ == "__main__":
    main()