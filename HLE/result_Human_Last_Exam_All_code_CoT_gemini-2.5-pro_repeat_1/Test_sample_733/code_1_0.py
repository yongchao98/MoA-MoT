def check_maximality(group_name, elements, op_func, s_candidate):
    """
    Checks if a set S is a maximal product-free set of size 2 in a group G.
    
    Args:
        group_name (str): The name of the group.
        elements (list): A list of the group's elements.
        op_func (function): The group operation.
        s_candidate (list): The candidate maximal product-free set.

    Returns:
        bool: True if S is a maximal product-free set, False otherwise.
    """
    # 1. Verify that S is product-free
    for x in s_candidate:
        for y in s_candidate:
            if op_func(x, y) in s_candidate:
                # This should not happen for our well-chosen sets
                return False

    # 2. Verify maximality: for any g not in S, S U {g} is not product-free
    for g in elements:
        if g not in s_candidate:
            s_extended = s_candidate + [g]
            
            # Check if s_extended is product-free. If it is, then s_candidate is not maximal.
            is_extension_product_free = True
            for x in s_extended:
                for y in s_extended:
                    if op_func(x, y) in s_extended:
                        is_extension_product_free = False
                        break
                if not is_extension_product_free:
                    break
            
            if is_extension_product_free:
                # Found a g that can be added, so s_candidate is not maximal
                return False

    # If all checks pass, the set is maximal product-free
    return True

def solve():
    """
    Finds and verifies the number of finite groups with a maximal product-free set of size 2.
    """
    # Define the 5 candidate groups and their respective sets
    
    # Group C4: {0, 1, 2, 3} with addition mod 4. Set S = {1, 3}
    c4 = {'name': 'C_4', 'elements': list(range(4)), 'op': lambda a, b: (a + b) % 4, 's': [1, 3]}

    # Group C2 x C2 (V4): {(0,0), (0,1), (1,0), (1,1)}. Set S = {(0,1), (1,0)}
    v4 = {'name': 'C_2 x C_2', 'elements': [(0,0), (0,1), (1,0), (1,1)], 
          'op': lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2), 's': [(0,1), (1,0)]}
          
    # Group C5: {0, 1, 2, 3, 4}. Set S = {2, 3}
    c5 = {'name': 'C_5', 'elements': list(range(5)), 'op': lambda a, b: (a + b) % 5, 's': [2, 3]}

    # Group C6: {0, 1, 2, 3, 4, 5}. Set S = {1, 4}
    c6 = {'name': 'C_6', 'elements': list(range(6)), 'op': lambda a, b: (a + b) % 6, 's': [1, 4]}

    # Group S3: {e, (12), (13), (23), (123), (132)} -> {0, 1, 2, 3, 4, 5}
    # Set S = {(12), (123)} -> {1, 4}
    s3_mult_table = [
        [0, 1, 2, 3, 4, 5], [1, 0, 4, 5, 2, 3], [2, 5, 0, 4, 3, 1],
        [3, 4, 5, 0, 1, 2], [4, 3, 1, 2, 5, 0], [5, 2, 3, 1, 0, 4]
    ]
    s3 = {'name': 'S_3', 'elements': list(range(6)), 
          'op': lambda a, b: s3_mult_table[a][b], 's': [1, 4]}
          
    groups_to_check = [c4, v4, c5, c6, s3]
    
    found_groups = []
    for group in groups_to_check:
        if check_maximality(group['name'], group['elements'], group['op'], group['s']):
            found_groups.append(group['name'])

    count = len(found_groups)
    equation_parts = ["1"] * count
    
    print(f"Found {count} groups with a maximal product-free set of size 2.")
    print(f"The groups are: {', '.join(found_groups)}")
    if count > 0:
        print("Calculation: " + " + ".join(equation_parts) + f" = {count}")

solve()