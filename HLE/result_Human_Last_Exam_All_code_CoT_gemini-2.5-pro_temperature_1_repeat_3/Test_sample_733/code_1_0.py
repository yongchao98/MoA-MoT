import itertools

def get_group_op(mult_table, elements):
    """Creates a group operation function from a Cayley (multiplication) table."""
    element_map = {el: i for i, el in enumerate(elements)}
    def op(a, b):
        i = element_map[a]
        j = element_map[b]
        return mult_table[i][j]
    return op

def is_product_free(s, op):
    """Checks if a subset S is product-free."""
    s_set = set(s)
    for x in s:
        for y in s:
            if op(x, y) in s_set:
                return False
    return True

def is_maximal(s, op, elements):
    """Checks if a product-free set S is maximal by inclusion."""
    s_set = set(s)
    # Get all elements of the group not in the set s
    g_minus_s = [el for el in elements if el not in s_set]
    
    # For s to be maximal, adding any element g from outside s should make it NOT product-free
    for g in g_minus_s:
        s_union_g = s + [g]
        # is_product_free is True if NO product falls back into the set.
        # If we find even one g for which s_union_g is still product-free, then s is not maximal.
        if is_product_free(s_union_g, op):
            return False
            
    # If all extensions are not product-free, then s is maximal
    return True

def find_groups_with_property():
    """
    Finds which non-isomorphic finite groups from a candidate list
    contain a maximal by inclusion product-free set of size 2.
    """
    found_groups = []
    
    # Define candidate groups structures (elements and Cayley tables/operations)
    groups = {
        "Z_3 (Cyclic group of order 3)": {
            "elements": [0, 1, 2],
            "op": lambda a, b: (a + b) % 3
        },
        "Z_4 (Cyclic group of order 4)": {
            "elements": [0, 1, 2, 3],
            "op": lambda a, b: (a + b) % 4
        },
        "V_4 (Klein four-group, Z_2 x Z_2)": {
            "elements": ['e', 'a', 'b', 'c'],
            "op": get_group_op([
                ['e', 'a', 'b', 'c'],
                ['a', 'e', 'c', 'b'],
                ['b', 'c', 'e', 'a'],
                ['c', 'b', 'a', 'e']
            ], ['e', 'a', 'b', 'c'])
        },
        "Z_5 (Cyclic group of order 5)": {
            "elements": [0, 1, 2, 3, 4],
            "op": lambda a, b: (a + b) % 5
        },
        "Z_6 (Cyclic group of order 6)": {
            "elements": [0, 1, 2, 3, 4, 5],
            "op": lambda a, b: (a + b) % 6
        },
        "S_3 (Symmetric group of order 6)": {
            "elements": ['e', 'r', 'r2', 's', 'sr', 'sr2'],
            "op": get_group_op([
                ['e', 'r', 'r2', 's', 'sr', 'sr2'], # e
                ['r', 'r2', 'e', 'sr', 'sr2', 's'], # r
                ['r2', 'e', 'r', 'sr2', 's', 'sr'], # r2
                ['s', 'sr2', 'sr', 'e', 'r2', 'r'], # s
                ['sr', 's', 'sr2', 'r', 'e', 'r2'], # sr
                ['sr2', 'sr', 's', 'r2', 'r', 'e']  # sr2
            ], ['e', 'r', 'r2', 's', 'sr', 'sr2'])
        },
        "Z_7 (Cyclic group of order 7)": {
            "elements": [0, 1, 2, 3, 4, 5, 6],
            "op": lambda a, b: (a + b) % 7
        }
    }

    print("Checking candidate groups...\n")
    for name, group_data in groups.items():
        elements = group_data["elements"]
        op = group_data["op"]
        
        has_property = False
        
        # Generate all subsets of size 2
        subsets_size_2 = itertools.combinations(elements, 2)
        
        for s_tuple in subsets_size_2:
            s = list(s_tuple)
            # Check if the set is product-free first
            if is_product_free(s, op):
                # If it is, check if it's maximal
                if is_maximal(s, op, elements):
                    has_property = True
                    break # Found one for this group, can move to the next group
        
        if has_property:
            found_groups.append(name)

    print("The finite groups that contain a maximal product-free set of size 2 are:")
    for name in sorted(found_groups):
        print(f"- {name}")
    
    count = len(found_groups)
    print(f"\nThe number of such finite groups is {count}.")


if __name__ == '__main__':
    find_groups_with_property()