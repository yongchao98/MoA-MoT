import itertools

def is_product_free(s, group_table, elements):
    """Checks if a set S is product-free."""
    s_indices = {elements.index(item) for item in s}
    for x_idx in s_indices:
        for y_idx in s_indices:
            product_idx = group_table[x_idx][y_idx]
            if elements[product_idx] in s:
                return False
    return True

def has_maximal_product_free_set_of_size_2(group_name, group_table):
    """
    Checks if a group has a maximal by inclusion product-free set of size 2.
    """
    n = len(group_table)
    elements = list(range(n))

    # Iterate over all subsets of size 2
    for s_tuple in itertools.combinations(elements, 2):
        s = set(s_tuple)

        # 1. Check if S is product-free
        if not is_product_free(s, group_table, elements):
            continue

        # 2. Check if S is maximal
        is_maximal = True
        g_minus_s = [g for g in elements if g not in s]
        
        for g in g_minus_s:
            s_prime = s.union({g})
            # If S_prime is product-free, then S is not maximal
            if is_product_free(s_prime, group_table, elements):
                is_maximal = False
                break
        
        if is_maximal:
            # Found one such set, so this group qualifies.
            # Example set found: S={a,b} where a,b are the elements at these indices.
            # print(f"Found maximal product-free set {s} in {group_name}")
            return True
            
    return False

def main():
    """
    Defines small groups by their Cayley tables and checks them.
    """
    groups = {
        # Trivial group
        "C1": [[0]],
        # Order 2
        "C2": [[0, 1], [1, 0]],
        # Order 3
        "C3": [[0, 1, 2], [1, 2, 0], [2, 0, 1]],
        # Order 4
        "C4": [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]],
        "V4 (C2xC2)": [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]],
        # Order 5
        "C5": [[0, 1, 2, 3, 4], [1, 2, 3, 4, 0], [2, 3, 4, 0, 1], [3, 4, 0, 1, 2], [4, 0, 1, 2, 3]],
        # Order 6
        "C6": [[(i + j) % 6 for j in range(6)] for i in range(6)],
        "S3": [[0, 1, 2, 3, 4, 5], 
               [1, 2, 0, 4, 5, 3], 
               [2, 0, 1, 5, 3, 4], 
               [3, 5, 4, 0, 2, 1], 
               [4, 3, 5, 1, 0, 2], 
               [5, 4, 3, 2, 1, 0]],
        # Order 7
        "C7": [[(i + j) % 7 for j in range(7)] for i in range(7)],
        # Order 8
        "Q8": [[0, 1, 2, 3, 4, 5, 6, 7],
               [1, 0, 3, 2, 5, 4, 7, 6],
               [2, 3, 1, 0, 6, 7, 5, 4],
               [3, 2, 0, 1, 7, 6, 4, 5],
               [4, 5, 7, 6, 1, 0, 2, 3],
               [5, 4, 6, 7, 0, 1, 3, 2],
               [6, 7, 4, 5, 3, 2, 1, 0],
               [7, 6, 5, 4, 2, 3, 0, 1]]
    }

    count = 0
    solution_groups = []
    for name, table in groups.items():
        if has_maximal_product_free_set_of_size_2(name, table):
            solution_groups.append(name)
            count += 1
            
    print("The finite groups found to have a maximal product-free set of size 2 are:")
    for name in solution_groups:
        print(f"- {name}")
    print(f"\nIn total, there are {count} such groups among those tested.")


if __name__ == "__main__":
    main()
