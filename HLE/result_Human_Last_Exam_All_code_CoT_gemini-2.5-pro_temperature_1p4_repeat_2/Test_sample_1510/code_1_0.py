def count_wallpaper_groups_by_point_group_order():
    """
    This script calculates the number of 2D wallpaper groups
    that have a point group of a specific order.
    """
    # Data of the 17 wallpaper groups with their point group and its order.
    wallpaper_groups_data = {
        'p1':   {'point_group': '1',   'order': 1},
        'p2':   {'point_group': '2',   'order': 2},
        'pm':   {'point_group': 'm',   'order': 2},
        'pg':   {'point_group': 'm',   'order': 2},
        'cm':   {'point_group': 'm',   'order': 2},
        'p2mm': {'point_group': '2mm', 'order': 4},
        'p2mg': {'point_group': '2mm', 'order': 4},
        'p2gg': {'point_group': '2mm', 'order': 4},
        'c2mm': {'point_group': '2mm', 'order': 4},
        'p4':   {'point_group': '4',   'order': 4},
        'p4mm': {'point_group': '4mm', 'order': 8},
        'p4gm': {'point_group': '4mm', 'order': 8},
        'p3':   {'point_group': '3',   'order': 3},
        'p3m1': {'point_group': '3m',  'order': 6},
        'p31m': {'point_group': '3m',  'order': 6},
        'p6':   {'point_group': '6',   'order': 6},
        'p6mm': {'point_group': '6mm', 'order': 12}
    }

    target_order = 2
    
    print(f"Finding wallpaper groups where the point group has an order of {target_order}.")
    print("-" * 60)

    # Group wallpaper groups by their point group name
    groups_by_pg = {}
    for wg_name, data in wallpaper_groups_data.items():
        if data['order'] == target_order:
            pg_name = data['point_group']
            if pg_name not in groups_by_pg:
                groups_by_pg[pg_name] = []
            groups_by_pg[pg_name].append(wg_name)
    
    counts = []
    total_count = 0

    # Print the breakdown
    for pg_name, wg_list in sorted(groups_by_pg.items()):
        count = len(wg_list)
        counts.append(str(count))
        total_count += count
        print(f"Point group '{pg_name}' has order {target_order}.")
        print(f"  - It corresponds to {count} wallpaper group(s): {', '.join(wg_list)}")
        
    # Print the final equation and total
    equation_str = " + ".join(counts)
    
    print("-" * 60)
    print(f"The calculation for the total is: {equation_str} = {total_count}")
    print(f"Therefore, there are a total of {total_count} wallpaper groups with a point group of order {target_order}.")

# Execute the function
count_wallpaper_groups_by_point_group_order()