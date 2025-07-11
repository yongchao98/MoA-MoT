import itertools

def is_product_free(table, subset):
    """Checks if a subset is product-free for a group given by its multiplication table."""
    n = len(table)
    # The identity element '0' cannot be in a product-free set.
    if 0 in subset:
        return False
        
    for s1 in subset:
        for s2 in subset:
            product = table[s1][s2]
            if product in subset:
                return False
    return True

def has_maximal_product_free_set_of_size_2(name, table):
    """
    Checks if a group has a maximal by inclusion product-free set of size 2.
    The identity element is assumed to be 0.
    """
    n = len(table)
    elements = range(1, n)  # Exclude identity element 0

    # Iterate through all 2-element subsets
    for subset_tuple in itertools.combinations(elements, 2):
        subset = set(subset_tuple)

        # 1. Check if the subset is product-free
        if not is_product_free(table, subset):
            continue

        # 2. If it is product-free, check if it's maximal
        is_maximal = True
        other_elements = [g for g in range(n) if g not in subset]
        
        for g in other_elements:
            # Check if S union {g} is product-free.
            # If we find a 'g' for which S union {g} is still product-free,
            # then S is not maximal.
            new_set = subset.union({g})
            if is_product_free(table, new_set):
                is_maximal = False
                break
        
        if is_maximal:
            # Found one such set, so the group qualifies.
            return True
            
    return False

def solve():
    """
    Finds all finite groups up to order 10 with a maximal product-free set of size 2.
    """
    # Multiplication tables for non-isomorphic groups up to order 10.
    # e = 0
    groups = {
        "C1": [[0]],
        "C2": [[0, 1], [1, 0]],
        "C3": [[0, 1, 2], [1, 2, 0], [2, 0, 1]],
        "C4": [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]],
        "V4": [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]],
        "C5": [[0, 1, 2, 3, 4], [1, 2, 3, 4, 0], [2, 3, 4, 0, 1], [3, 4, 0, 1, 2], [4, 0, 1, 2, 3]],
        "C6": [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0], [2, 3, 4, 5, 0, 1], [3, 4, 5, 0, 1, 2], [4, 5, 0, 1, 2, 3], [5, 0, 1, 2, 3, 4]],
        "S3": [[0, 1, 2, 3, 4, 5], [1, 2, 0, 4, 5, 3], [2, 0, 1, 5, 3, 4], [3, 5, 4, 0, 2, 1], [4, 3, 5, 1, 0, 2], [5, 4, 3, 2, 1, 0]],
        "C7": [[0,1,2,3,4,5,6],[1,2,3,4,5,6,0],[2,3,4,5,6,0,1],[3,4,5,6,0,1,2],[4,5,6,0,1,2,3],[5,6,0,1,2,3,4],[6,0,1,2,3,4,5]],
        "C8": [[0,1,2,3,4,5,6,7],[1,2,3,4,5,6,7,0],[2,3,4,5,6,7,0,1],[3,4,5,6,7,0,1,2],[4,5,6,7,0,1,2,3],[5,6,7,0,1,2,3,4],[6,7,0,1,2,3,4,5],[7,0,1,2,3,4,5,6]],
        "C4xC2":[[0,1,2,3,4,5,6,7],[1,2,3,0,5,6,7,4],[2,3,0,1,6,7,4,5],[3,0,1,2,7,4,5,6],[4,5,6,7,0,1,2,3],[5,6,7,4,1,2,3,0],[6,7,4,5,2,3,0,1],[7,4,5,6,3,0,1,2]],
        "C2x3": [[0,1,2,3,4,5,6,7],[1,0,3,2,5,4,7,6],[2,3,0,1,6,7,4,5],[3,2,1,0,7,6,5,4],[4,5,6,7,0,1,2,3],[5,4,7,6,1,0,3,2],[6,7,4,5,2,3,0,1],[7,6,5,4,3,2,1,0]],
        "D8": [[0,1,2,3,4,5,6,7],[1,2,3,0,5,6,7,4],[2,3,0,1,6,7,4,5],[3,0,1,2,7,4,5,6],[4,7,6,5,0,3,2,1],[5,4,7,6,1,0,3,2],[6,5,4,7,2,1,0,3],[7,6,5,4,3,2,1,0]],
        "Q8": [[0,1,2,3,4,5,6,7],[1,2,3,0,5,6,7,4],[2,3,0,1,6,7,4,5],[3,0,1,2,7,4,5,6],[4,7,6,5,1,3,2,0],[5,4,7,6,2,0,3,1],[6,5,4,7,3,1,0,2],[7,6,5,4,0,2,1,3]],
        "C9": [[(i+j)%9 for j in range(9)] for i in range(9)],
        "C3xC3": [[( (i//3+j//3)%3 )*3 + (i%3+j%3)%3 for j in range(9)] for i in range(9)],
        "C10": [[(i+j)%10 for j in range(10)] for i in range(10)],
        "D10": [[j if i < 5 else (j+i)%5+5 if j < 5 else 5-(j-i)%5 if (j-i)%5>0 else 0 for j in range(10)] for i in range(10)],
    }
    
    # In D10 table, e=0,r=1,r2=2,r3=3,r4=4, s=5,sr=6,sr2=7,sr3=8,sr4=9
    # a,b < 5: a+b %5. a>=5,b<5: (b+a-5)%5+5. a<5,b>=5: (a+b-5)%5+5. a,b>=5: (a-b)%5
    # My D10 table is slightly off, lets build it right.
    # r^i * r^j = r^(i+j), s*r^i=r^-i*s, r^i*s*r^j*s = r^(i-j)
    D10_table = [[0]*10 for _ in range(10)]
    for i in range(10):
        for j in range(10):
            is_i_s = i >= 5
            is_j_s = j >= 5
            i_r = i % 5
            j_r = j % 5
            if not is_i_s and not is_j_s: # r*r
                D10_table[i][j] = (i_r + j_r) % 5
            elif not is_i_s and is_j_s: # r*s
                D10_table[i][j] = ((i_r + j_r) % 5) + 5
            elif is_i_s and not is_j_s: # s*r
                D10_table[i][j] = ((-i_r + j_r) % 5) + 5
            else: # s*s
                D10_table[i][j] = ((-i_r + j_r) % 5)
    groups["D10"] = D10_table


    found_groups = []
    for name, table in groups.items():
        if len(table) < 3: # Need at least 3 elements for a size 2 subset + identity
            continue
        if has_maximal_product_free_set_of_size_2(name, table):
            found_groups.append(name)
            
    print("Groups containing a maximal by inclusion product-free set of size 2:")
    for name in sorted(found_groups):
        print(f"- {name}")
    
    count = len(found_groups)
    print(f"\nTotal number of such groups is: {count}")

solve()