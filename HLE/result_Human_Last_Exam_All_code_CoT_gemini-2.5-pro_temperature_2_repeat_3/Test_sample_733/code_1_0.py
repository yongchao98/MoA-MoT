import itertools

def is_product_free(s, group_elements, mul_table):
    """Checks if a set s is product-free."""
    s_indices = {group_elements.index(e) for e in s}
    for e1 in s:
        for e2 in s:
            idx1 = group_elements.index(e1)
            idx2 = group_elements.index(e2)
            product_idx = mul_table[idx1][idx2]
            if product_idx in s_indices:
                return False
    return True

def is_maximal(s, group_elements, mul_table):
    """Checks if a product-free set s is maximal."""
    s_with_g = set(s)
    other_elements = [e for e in group_elements if e not in s]
    for g in other_elements:
        s_with_g.add(g)
        # Check if S u {g} is product-free. If it is, then S is not maximal.
        if is_product_free(s_with_g, group_elements, mul_table):
            return False
        s_with_g.remove(g)
    return True

def find_groups_with_maximal_pf_set_of_size_2():
    """
    Tests a list of small finite groups and counts how many have
    a maximal by inclusion product-free set of size 2.
    """
    groups_to_test = {
        'C2': {
            'elements': ['e', 'a'],
            'table': [[0, 1], [1, 0]]
        },
        'C3': {
            'elements': ['e', 'a', 'a^2'],
            'table': [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
        },
        'C4': {
            'elements': ['e', 'a', 'a^2', 'a^3'],
            'table': [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]
        },
        'V4': { # Klein four-group (C2 x C2)
            'elements': ['e', 'a', 'b', 'c'],
            'table': [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
        },
        'C5': {
            'elements': ['e', 'a', 'a^2', 'a^3', 'a^4'],
            'table': [[0, 1, 2, 3, 4], [1, 2, 3, 4, 0], [2, 3, 4, 0, 1],
                      [3, 4, 0, 1, 2], [4, 0, 1, 2, 3]]
        },
        'C6': {
             'elements': ['e', 'a', 'a^2', 'a^3', 'a^4', 'a^5'],
             'table': [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0], [2, 3, 4, 5, 0, 1],
                       [3, 4, 5, 0, 1, 2], [4, 5, 0, 1, 2, 3], [5, 0, 1, 2, 3, 4]]
        },
        'S3': { # Symmetric group on 3 elements
            'elements': ['e', 'r', 'r^2', 's', 'sr', 'sr^2'],
            # e, r, r2, s, sr, sr2  (indices 0-5)
            # r^3=e, s^2=e, rs=sr^2
            'table': [[0, 1, 2, 3, 4, 5], # e
                      [1, 2, 0, 5, 3, 4], # r
                      [2, 0, 1, 4, 5, 3], # r^2
                      [3, 4, 5, 0, 1, 2], # s
                      [4, 5, 3, 2, 0, 1], # sr
                      [5, 3, 4, 1, 2, 0]] # sr^2
        },
         'C7': {
            'elements': ['e', 'a', 'a^2', 'a^3', 'a^4', 'a^5', 'a^6'],
            'table': [[0, 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6, 0], [2, 3, 4, 5, 6, 0, 1],
                      [3, 4, 5, 6, 0, 1, 2], [4, 5, 6, 0, 1, 2, 3], [5, 6, 0, 1, 2, 3, 4],
                      [6, 0, 1, 2, 3, 4, 5]]
        },
    }

    found_groups_count = 0
    print("Checking for finite groups with maximal by inclusion product-free sets of size 2...")
    print("-" * 20)
    for name, group_data in groups_to_test.items():
        elements = group_data['elements']
        mul_table = group_data['table']
        non_identity_elements = [e for e in elements if e != 'e']
        
        if len(non_identity_elements) < 2:
            continue
            
        found_in_group = False
        # Generate all subsets of size 2 from non-identity elements
        for s_tuple in itertools.combinations(non_identity_elements, 2):
            s = set(s_tuple)
            if is_product_free(s, elements, mul_table):
                if is_maximal(s, elements, mul_table):
                    found_groups_count += 1
                    print(f"Found a qualifying set in group {name}. Set: {s}")
                    found_in_group = True
                    break # Move to the next group once one set is found
        if not found_in_group:
             print(f"Group {name} does not have such a set.")
        print("-" * 20)

    print(f"\nTotal number of such groups found: {found_groups_count}")
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # Since this is a counting problem, there is no literal equation. I will show the sum.
    components = ['1'] * found_groups_count
    equation_str = " + ".join(components)
    if not equation_str:
      equation_str = "0"
    print(f"Final calculation: {equation_str} = {found_groups_count}")


if __name__ == '__main__':
    find_groups_with_maximal_pf_set_of_size_2()