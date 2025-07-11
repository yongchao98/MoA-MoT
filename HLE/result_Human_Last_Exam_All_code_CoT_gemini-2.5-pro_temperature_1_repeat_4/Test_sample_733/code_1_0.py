import itertools

def is_product_free(s, n, cayley_table):
    """Checks if a set s is product-free."""
    elements = list(s)
    for s1 in elements:
        for s2 in elements:
            product = cayley_table[s1][s2]
            if product in s:
                return False
    return True

def is_maximal_product_free(s, n, cayley_table):
    """Checks if a set s is a maximal product-free set."""
    if not is_product_free(s, n, cayley_table):
        return False
    
    all_elements = set(range(n))
    other_elements = all_elements - s
    
    for g in other_elements:
        s_prime = s.union({g})
        if is_product_free(s_prime, n, cayley_table):
            # Found a larger product-free set, so s is not maximal
            return False
            
    return True

def has_maximal_product_free_set_of_size_2(n, cayley_table):
    """Checks if a group has a maximal product-free set of size 2."""
    elements = range(n)
    for s_tuple in itertools.combinations(elements, 2):
        s = set(s_tuple)
        if is_maximal_product_free(s, n, cayley_table):
            return True
    return False

def get_cayley_table(group_name):
    """Returns the order and Cayley table for a given group."""
    if group_name == "C2":
        n = 2
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "C3":
        n = 3
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "C4":
        n = 4
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "C2xC2":
        n = 4
        # Elements 0, 1, 2, 3 correspond to (0,0), (1,0), (0,1), (1,1)
        # Operation is component-wise XOR
        table = [[i ^ j for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "C5":
        n = 5
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "C6":
        n = 6
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    if group_name == "S3":
        n = 6
        # e, (123), (132), (12), (13), (23) -> 0, 1, 2, 3, 4, 5
        table = [
            [0, 1, 2, 3, 4, 5],
            [1, 2, 0, 4, 5, 3],
            [2, 0, 1, 5, 3, 4],
            [3, 5, 4, 0, 2, 1],
            [4, 3, 5, 1, 0, 2],
            [5, 4, 3, 2, 1, 0]
        ]
        return n, table
    if group_name == "C7":
        n = 7
        table = [[(i + j) % n for j in range(n)] for i in range(n)]
        return n, table
    return 0, []

if __name__ == '__main__':
    groups_to_check = ["C2", "C3", "C4", "C2xC2", "C5", "C6", "S3", "C7"]
    
    found_groups = []
    
    for name in groups_to_check:
        n, table = get_cayley_table(name)
        if n > 0 and has_maximal_product_free_set_of_size_2(n, table):
            found_groups.append(name)
            
    print("The finite groups containing maximal by inclusion product-free sets of size 2 are:")
    for name in found_groups:
        print(f"- {name}")
    
    count = len(found_groups)
    
    equation_parts = ["1"] * count
    equation_str = " + ".join(equation_parts)
    
    print("\nFinal calculation:")
    if count > 0:
        print(f"{equation_str} = {count}")
    else:
        print("0")
