def count_wallpaper_groups_by_point_group_order():
    """
    This script calculates the number of wallpaper groups that have a point group of order 2.
    """

    # A mapping of the 17 wallpaper groups to their corresponding crystallographic point group.
    wallpaper_to_point_group = {
        'p1': '1',
        'p2': '2',
        'pm': 'm',
        'pg': 'm',
        'cm': 'm',
        'p2mm': '2mm',
        'p2mg': '2mm',
        'p2gg': '2mm',
        'c2mm': '2mm',
        'p3': '3',
        'p3m1': '3m',
        'p31m': '3m',
        'p4': '4',
        'p4mm': '4mm',
        'p4gm': '4mm',
        'p6': '6',
        'p6mm': '6mm',
    }

    # The order (number of symmetry operations) for each point group.
    point_group_orders = {
        '1': 1,
        '2': 2,    # Point group of order 2
        'm': 2,    # Point group of order 2
        '2mm': 4,
        '3': 3,
        '3m': 6,
        '4': 4,
        '4mm': 8,
        '6': 6,
        '6mm': 12,
    }

    target_order = 2
    groups_by_pg = {}

    # Iterate through all wallpaper groups and categorize them by their point group
    # if the point group has the target order.
    for wg_name, pg_name in wallpaper_to_point_group.items():
        if point_group_orders.get(pg_name) == target_order:
            if pg_name not in groups_by_pg:
                groups_by_pg[pg_name] = []
            groups_by_pg[pg_name].append(wg_name)
    
    print("Finding wallpaper groups with a point group of order 2...")
    
    counts = []
    total_count = 0

    # Print the details for each relevant point group
    for pg_name, wg_list in sorted(groups_by_pg.items()):
        count = len(wg_list)
        print(f"\nPoint group '{pg_name}' has order {target_order}.")
        print(f"Associated wallpaper groups: {', '.join(sorted(wg_list))}")
        print(f"Count: {count}")
        counts.append(str(count))
        total_count += count
    
    # Print the final equation as requested
    final_equation = " + ".join(counts)
    print(f"\nThe total number is the sum of these counts.")
    print(f"Total = {final_equation} = {total_count}")

# Execute the function
count_wallpaper_groups_by_point_group_order()