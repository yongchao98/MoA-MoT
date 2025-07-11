from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.named_groups import (
    SymmetricGroup, DihedralGroup, AlternatingGroup, CyclicGroup
)
from sympy.combinatorics.group_constructs import direct_product
from itertools import combinations

def get_group_info(G):
    """
    Generates elements, an element-to-index map, a Cayley table,
    and the index of the identity element for a given SymPy group.
    """
    # Use generate_dimino() for potentially faster generation for small groups
    elements = list(G.generate_dimino())
    n = len(elements)
    element_map = {el: i for i, el in enumerate(elements)}
    
    table = [[0] * n for _ in range(n)]
    for i, p1 in enumerate(elements):
        for j, p2 in enumerate(elements):
            prod = p1 * p2
            table[i][j] = element_map[prod]
            
    identity_index = element_map[G.identity]
    return elements, element_map, table, identity_index

def is_product_free(s_indices, cayley_table):
    """
    Checks if a set of indices corresponds to a product-free set.
    """
    for i in s_indices:
        for j in s_indices:
            product_index = cayley_table[i][j]
            if product_index in s_indices:
                return False
    return True

def has_maximal_product_free_set_of_size_3(G):
    """
    Checks if a group G has a maximal product-free set of size 3.
    """
    # Group must have at least 4 elements: identity + 3 for the set.
    if G.order() < 4:
        return False
        
    elements, _, cayley_table, identity_index = get_group_info(G)
    
    element_indices = list(range(len(elements)))
    candidate_indices = [i for i in element_indices if i != identity_index]
    
    # Must have at least 3 non-identity elements.
    if len(candidate_indices) < 3:
        return False

    # Iterate through all 3-element subsets of non-identity elements.
    for s_indices in combinations(candidate_indices, 3):
        if is_product_free(s_indices, cayley_table):
            # The set S is product-free. Now check for maximality.
            is_maximal = True
            other_indices = [i for i in element_indices if i not in s_indices]
            for g_index in other_indices:
                # If S can be extended by g, it's not maximal.
                T_indices = list(s_indices) + [g_index]
                if is_product_free(T_indices, cayley_table):
                    is_maximal = False
                    break 
            if is_maximal:
                # Found a maximal product-free set of size 3.
                return True
                
    return False

def find_groups_and_count():
    """
    Tests a list of finite groups and prints the total count.
    """
    # List of non-isomorphic groups up to order 12, plus a few from literature.
    # Note: S3xZ2 is isomorphic to D12
    groups_to_check = {
        "Z4": CyclicGroup(4),
        "V4": direct_product(CyclicGroup(2), CyclicGroup(2)),
        "Z5": CyclicGroup(5),
        "Z6": CyclicGroup(6),
        "S3": SymmetricGroup(3),
        "Z7": CyclicGroup(7),
        "Z8": CyclicGroup(8),
        "Z4xZ2": direct_product(CyclicGroup(4), CyclicGroup(2)),
        "Z2xZ2xZ2": direct_product(CyclicGroup(2), CyclicGroup(2), CyclicGroup(2)),
        "D8": DihedralGroup(4),  # Order 8
        "Z9": CyclicGroup(9),
        "Z3xZ3": direct_product(CyclicGroup(3), CyclicGroup(3)),
        "Z10": CyclicGroup(10),
        "D10": DihedralGroup(5), # Order 10
        "A4": AlternatingGroup(4), # Order 12
        "D12": DihedralGroup(6), # Order 12
    }
    
    # Q8 must be constructed manually as it is not a named group in SymPy.
    i = Permutation(0, 1, 2, 3)(4, 5, 6, 7)
    j = Permutation(0, 4, 2, 6)(1, 7, 3, 5)
    Q8 = i.group_closure([i, j])
    groups_to_check["Q8"] = Q8

    found_groups = []
    for name, G in groups_to_check.items():
        if has_maximal_product_free_set_of_size_3(G):
            found_groups.append(name)
            
    print("Based on a computational search using standard definitions, the following groups contain a maximal by inclusion product-free set of size 3:")
    for name in found_groups:
        print(f"- {name}")
    print("\nTotal number of such groups found:")
    print(len(found_groups))

if __name__ == '__main__':
    find_groups_and_count()