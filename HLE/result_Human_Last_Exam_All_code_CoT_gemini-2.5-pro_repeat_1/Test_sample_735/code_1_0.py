import itertools

def is_product_free(s, group_table):
    """Checks if a subset s is product-free."""
    s_indices = {elem for elem in s}
    for i in s:
        for j in s:
            product = group_table[i][j]
            if product in s_indices:
                return False
    return True

def find_maximal_product_free_sets_of_size_3(group_name, elements, group_table):
    """
    Finds if a group contains a maximal product-free set of size 3.
    """
    n = len(elements)
    if n < 3:
        return False
        
    element_map = {name: i for i, name in enumerate(elements)}
    indices = list(range(n))

    # Iterate over all subsets of size 3
    for s_indices in itertools.combinations(indices, 3):
        if is_product_free(s_indices, group_table):
            # Check for maximality
            is_maximal = True
            other_indices = [i for i in indices if i not in s_indices]
            
            for g_idx in other_indices:
                # Test if S U {g} is product-free
                s_prime = list(s_indices) + [g_idx]
                if is_product_free(s_prime, group_table):
                    # Found a larger product-free set, so s is not maximal
                    is_maximal = False
                    break
            
            if is_maximal:
                # Found one, so this group qualifies.
                return True
    return False

def main():
    """
    Main function to define small groups and check them.
    """
    # Cayley tables for small groups. Elements are represented by their indices 0, 1, 2, ...
    groups = {
        "Z_1": ([0], [[0]]),
        "Z_2": ([0, 1], [[0, 1], [1, 0]]),
        "Z_3": ([0, 1, 2], [[0, 1, 2], [1, 2, 0], [2, 0, 1]]),
        "Z_4": ([0, 1, 2, 3], [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]),
        "V4 (Z_2xZ_2)": ([0, 1, 2, 3], [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]),
        "Z_5": ([0, 1, 2, 3, 4], [[(i + j) % 5 for j in range(5)] for i in range(5)]),
        "Z_6": ([0, 1, 2, 3, 4, 5], [[(i + j) % 6 for j in range(6)] for i in range(6)]),
        "S_3": ([0, 1, 2, 3, 4, 5], [ # e, a, a^2, b, ab, a^2b
            [0, 1, 2, 3, 4, 5], [1, 2, 0, 4, 5, 3], [2, 0, 1, 5, 3, 4],
            [3, 5, 4, 0, 2, 1], [4, 3, 5, 1, 0, 2], [5, 4, 3, 2, 1, 0]
        ]),
        "Z_7": ([0,1,2,3,4,5,6], [[(i+j)%7 for j in range(7)] for i in range(7)]),
        "Z_8": ([0,1,2,3,4,5,6,7], [[(i+j)%8 for j in range(8)] for i in range(8)]),
        "Z_9": ([0,1,2,3,4,5,6,7,8], [[(i+j)%9 for j in range(9)] for i in range(9)]),
        "Z_3xZ_3": (
            [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)],
            [[( (i//3+j//3)%3, (i%3+j%3)%3 ) for j in range(9)] for i in range(9)]
        )
    }

    # The Cayley table for Z_3xZ_3 needs to be mapped back to indices
    z3z3_elements, z3z3_tuples = groups["Z_3xZ_3"]
    z3z3_map = {name: i for i, name in enumerate(z3z3_elements)}
    z3z3_table = [[z3z3_map[val] for val in row] for row in z3z3_tuples]
    groups["Z_3xZ_3"] = (list(range(9)), z3z3_table)
    
    found_groups = []
    for name, (elements, table) in groups.items():
        if find_maximal_product_free_sets_of_size_3(name, elements, table):
            found_groups.append(name)
            
    print("Found groups containing a maximal product-free set of size 3:")
    for name in found_groups:
        print(f"- {name}")
    
    count = len(found_groups)
    print(f"\nFinal Answer: The number of finite groups containing such a set is {count}.")


if __name__ == '__main__':
    main()