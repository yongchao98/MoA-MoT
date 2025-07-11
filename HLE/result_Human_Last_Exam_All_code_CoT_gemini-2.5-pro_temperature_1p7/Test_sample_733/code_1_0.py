import itertools

def is_product_free(s, cayley_table):
    """Checks if a set s is product-free given a cayley_table."""
    s_set = set(s)
    for x in s:
        for y in s:
            product = cayley_table[x][y]
            if product in s_set:
                return False
    return True

def find_maximal_product_free_set_of_size_2(group_name, cayley_table):
    """
    Finds if a group contains a maximal by inclusion product-free set of size 2.
    The group elements are represented as integers from 0 to n-1.
    """
    n = len(cayley_table)
    elements = list(range(n))
    
    # Iterate through all subsets of size 2
    for s in itertools.combinations(elements, 2):
        s_list = list(s)
        
        # 1. Check if the subset is product-free
        if is_product_free(s_list, cayley_table):
            
            # 2. If it is, check if it's maximal
            is_maximal = True
            other_elements = [el for el in elements if el not in s_list]
            
            for g in other_elements:
                extended_set = s_list + [g]
                # If we can extend the set and it's still product-free, then s is not maximal
                if is_product_free(extended_set, cayley_table):
                    is_maximal = False
                    break
            
            if is_maximal:
                # Found one, so this group qualifies.
                # We can print the set we found for verification.
                # print(f"  - Found maximal product-free set {s_list} in {group_name}")
                return True
                
    return False

def get_cayley_tables():
    """Returns a dictionary of group names to their Cayley tables."""
    groups = {}

    # C3 (Cyclic group of order 3)
    groups['C3'] = [[(i + j) % 3 for j in range(3)] for i in range(3)]
    
    # C4 (Cyclic group of order 4)
    groups['C4'] = [[(i + j) % 4 for j in range(4)] for i in range(4)]
    
    # C5 (Cyclic group of order 5)
    groups['C5'] = [[(i + j) % 5 for j in range(5)] for i in range(5)]
    
    # C2xC2 (Klein four-group)
    groups['C2xC2'] = [[i ^ j for j in range(4)] for i in range(4)]
    
    # D10 (Dihedral group of order 10)
    # Elements 0..4 are r^i, 5..9 are sr^i
    d10_table = [[0]*10 for _ in range(10)]
    for i in range(10):
        for j in range(10):
            if i < 5 and j < 5: d10_table[i][j] = (i + j) % 5      # r^i * r^j = r^(i+j)
            elif i < 5 and j >= 5: d10_table[i][j] = 5 + (j - 5 - i + 5) % 5 # r^i * sr^j = sr^(j-i)
            elif i >= 5 and j < 5: d10_table[i][j] = 5 + (i - 5 + j) % 5   # sr^i * r^j = sr^(i+j)
            else: d10_table[i][j] = (j - 5 - (i - 5) + 5) % 5 # sr^i * sr^j = r^(j-i)
    groups['D10'] = d10_table

    # Q12 (Dicyclic group of order 12)
    # Elements 0..5 are a^k, 6..11 are xa^k
    q12_table = [[0]*12 for _ in range(12)]
    for i in range(12):
        for j in range(12):
            if i < 6 and j < 6: q12_table[i][j] = (i + j) % 6             # a^i * a^j = a^(i+j)
            elif i < 6 and j >= 6: q12_table[i][j] = 6 + (j - 6 - i + 6) % 6 # a^i * xa^j = xa^(j-i)
            elif i >= 6 and j < 6: q12_table[i][j] = 6 + (i - 6 + j) % 6      # xa^i * a^j = xa^(i+j)
            else: q12_table[i][j] = (i - 6 + j - 6 + 3) % 6 # xa^i * xa^j = a^(i+j+3)
    groups['Q12'] = q12_table
    
    return groups

if __name__ == "__main__":
    cayley_tables = get_cayley_tables()
    found_groups_count = 0
    
    print("Checking for finite groups with maximal product-free sets of size 2...")
    
    for name, table in cayley_tables.items():
        if find_maximal_product_free_set_of_size_2(name, table):
            found_groups_count += 1
            print(f"Found qualifying group: {name}")

    print("\n---")
    print(f"The total number of such groups found is: {found_groups_count}")
