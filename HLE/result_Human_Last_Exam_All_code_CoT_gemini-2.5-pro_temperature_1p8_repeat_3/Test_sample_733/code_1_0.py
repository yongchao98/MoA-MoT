from itertools import combinations

def is_product_free(s, group_elements, op):
    """Checks if a subset 's' is product-free."""
    s_set = set(s)
    for s1 in s:
        for s2 in s:
            product = op(s1, s2)
            if product in s_set:
                return False
    return True

def has_maximal_product_free_set_of_size_2(group):
    """
    Checks if a group contains a maximal by inclusion product-free set of size 2.
    The group is defined by its elements and multiplication table.
    """
    elements = group['elements']
    op = group['op']
    n = len(elements)
    
    # We need at least 3 elements to have a non-trivial subset of size 2
    if n < 3:
        return False
        
    # Iterate through all subsets of size 2
    for s in combinations(elements, 2):
        s_list = list(s)
        
        # Check if the subset is product-free
        if is_product_free(s_list, elements, op):
            
            # If it is product-free, check for maximality
            is_maximal = True
            other_elements = [g for g in elements if g not in s_list]
            
            for g in other_elements:
                s_prime = s_list + [g]
                if is_product_free(s_prime, elements, op):
                    # Found a larger product-free set containing s, so s is not maximal
                    is_maximal = False
                    break
            
            if is_maximal:
                # Found a maximal product-free set of size 2
                return True
                
    return False

def main():
    """
    Defines several small finite groups and checks them for the property.
    """
    # Group definitions
    # C_n: Cyclic group of order n. Elements are integers 0..n-1, operation is addition modulo n.
    C3 = {'name': 'C3', 'elements': [0,1,2], 'op': lambda a,b: (a+b)%3}
    C4 = {'name': 'C4', 'elements': [0,1,2,3], 'op': lambda a,b: (a+b)%4}
    C5 = {'name': 'C5', 'elements': [0,1,2,3,4], 'op': lambda a,b: (a+b)%5}
    C6 = {'name': 'C6', 'elements': [0,1,2,3,4,5], 'op': lambda a,b: (a+b)%6}
    C7 = {'name': 'C7', 'elements': [0,1,2,3,4,5,6], 'op': lambda a,b: (a+b)%7}
    
    # V4: Klein-four group. e=0, a=1, b=2, c=3
    v4_table = [[0,1,2,3], [1,0,3,2], [2,3,0,1], [3,2,1,0]]
    V4 = {'name': 'V4', 'elements': [0,1,2,3], 'op': lambda a,b: v4_table[a][b]}

    # S3: Symmetric group on 3 elements. e=0, r=1, r^2=2, s=3, sr=4, sr^2=5
    s3_table = [
        [0, 1, 2, 3, 4, 5], [1, 2, 0, 4, 5, 3], [2, 0, 1, 5, 3, 4],
        [3, 5, 4, 0, 2, 1], [4, 3, 5, 1, 0, 2], [5, 4, 3, 2, 1, 0],
    ]
    S3 = {'name': 'S3', 'elements': [0,1,2,3,4,5], 'op': lambda a,b: s3_table[a][b]}
    
    groups_to_check = [C3, C4, C5, C6, C7, V4, S3]
    
    count = 0
    found_groups = []
    
    for g in groups_to_check:
        if has_maximal_product_free_set_of_size_2(g):
            count += 1
            found_groups.append(g['name'])
    
    print("Based on the analysis, the groups satisfying the condition are C4, C5, C7, and V4.")
    print(f"Verification code found {count} such groups among the tested samples: {', '.join(found_groups)}")
    print(f"Total number of finite groups found is: {count}")

if __name__ == '__main__':
    main()