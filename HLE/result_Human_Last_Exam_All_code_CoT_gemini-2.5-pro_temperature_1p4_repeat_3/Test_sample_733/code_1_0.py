from itertools import combinations

def generate_cayley_table(group_name):
    """Generates the Cayley (multiplication) table for a given group."""
    if group_name.startswith("Z"):
        # Cyclic group Z_n under addition modulo n
        try:
            n = int(group_name[1:])
            return [[(i + j) % n for j in range(n)] for i in range(n)]
        except (ValueError, IndexError):
            return None
            
    if group_name.startswith("D"):
        # Dihedral group D_{2n} of order 2n
        # Elements 0 to n-1 are rotations r^i
        # Elements n to 2n-1 are reflections sr^i
        try:
            n = int(group_name[1:]) // 2
            if 2 * n != int(group_name[1:]): return None
            order = 2 * n
            table = [[0] * order for _ in range(order)]
            # r^i * r^j = r^(i+j)
            for i in range(n):
                for j in range(n):
                    table[i][j] = (i + j) % n
            # r^i * sr^j = sr^(j-i)
            for i in range(n):
                for j in range(n):
                    table[i][n + j] = n + ((j - i) % n)
            # sr^i * r^j = sr^(i+j)
            for i in range(n):
                for j in range(n):
                    table[n + i][j] = n + ((i + j) % n)
            # sr^i * sr^j = r^(j-i)
            for i in range(n):
                for j in range(n):
                    table[n + i][n + j] = (j - i) % n
            return table
        except (ValueError, IndexError):
            return None

    return None

def is_product_free(S, table):
    """Checks if a set S is product-free given a Cayley table."""
    for x in S:
        for y in S:
            product = table[x][y]
            if product in S:
                return False
    return True

def has_maximal_product_free_set_of_size_2(table):
    """
    Checks if a group (represented by its table) has a maximal
    product-free set of size 2.
    """
    if not table:
        return False
        
    n = len(table)
    elements = list(range(n))

    # Iterate through all subsets of size 2
    for s_tuple in combinations(elements, 2):
        S = set(s_tuple)

        # 1. Check if S is product-free
        if not is_product_free(S, table):
            continue

        # 2. Check if S is maximal
        is_maximal = True
        other_elements = [g for g in elements if g not in S]
        for g in other_elements:
            S_prime = S.union({g})
            # If S_prime is product-free, then S is not maximal
            if is_product_free(S_prime, table):
                is_maximal = False
                break
        
        if is_maximal:
            # Found one for this group, so the group qualifies.
            return True
            
    return False

def main():
    """
    Main function to find the number of specified groups.
    We test the groups classified by G.L. Walls (2007) and some others as counter-examples.
    """
    groups_to_test = [
        "Z3", "Z4", "Z5", "Z6", "Z7",  # Cyclic groups
        "D6", "D8", "D10"             # Dihedral groups (D6=S3)
    ]
    
    found_groups = []
    
    print("Checking various finite groups...")
    for group_name in groups_to_test:
        table = generate_cayley_table(group_name)
        if has_maximal_product_free_set_of_size_2(table):
            found_groups.append(group_name)
            
    print("\nThe finite groups found to contain a maximal by inclusion product-free set of size 2 are:")
    for name in found_groups:
        print(f"- {name}")
    
    print(f"\nThe total number of such non-isomorphic groups found is:")
    print(len(found_groups))

if __name__ == "__main__":
    main()