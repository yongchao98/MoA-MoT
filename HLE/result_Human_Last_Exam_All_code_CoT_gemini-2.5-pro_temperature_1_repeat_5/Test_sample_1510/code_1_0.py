import collections

def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups that have a point group of a specific order.
    """
    # A dictionary mapping each of the 17 wallpaper groups to the order of its point group.
    # The point group (or crystal class) determines the rotational and reflectional symmetry at a point.
    # C_n has order n. D_n has order 2n.
    # C1: 1, C2: 2, C3: 3, C4: 4, C6: 6
    # D1: 2, D2: 4, D3: 6, D4: 8, D6: 12
    wallpaper_groups_data = {
        'p1': {'point_group': 'C1', 'order': 1},
        'p2': {'point_group': 'C2', 'order': 2},
        'pm': {'point_group': 'D1', 'order': 2},
        'pg': {'point_group': 'D1', 'order': 2},
        'cm': {'point_group': 'D1', 'order': 2},
        'p2mm': {'point_group': 'D2', 'order': 4},
        'p2mg': {'point_group': 'D2', 'order': 4},
        'p2gg': {'point_group': 'D2', 'order': 4},
        'c2mm': {'point_group': 'D2', 'order': 4},
        'p3': {'point_group': 'C3', 'order': 3},
        'p3m1': {'point_group': 'D3', 'order': 6},
        'p31m': {'point_group': 'D3', 'order': 6},
        'p4': {'point_group': 'C4', 'order': 4},
        'p4mm': {'point_group': 'D4', 'order': 8},
        'p4gm': {'point_group': 'D4', 'order': 8},
        'p6': {'point_group': 'C6', 'order': 6},
        'p6mm': {'point_group': 'D6', 'order': 12},
    }

    target_order = 2
    found_groups = []

    print(f"Finding wallpaper groups with a point group of order {target_order}:")
    for name, data in wallpaper_groups_data.items():
        if data['order'] == target_order:
            found_groups.append(name)
            print(f"- Found group '{name}' (Point Group: {data['point_group']})")

    count = len(found_groups)
    
    # Create the equation string as requested
    equation_parts = ['1'] * count
    equation_str = " + ".join(equation_parts)

    print(f"\nIn total, there are {count} such groups.")
    print(f"The final calculation is: {equation_str} = {count}")

count_wallpaper_groups_by_point_group_order()