def solve_wallpaper_groups():
    """
    This function calculates and explains the number of wallpaper groups
    with a point group of order 2.
    """
    # A mapping of each wallpaper group to its corresponding point group.
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
        'p4': '4',
        'p4mm': '4mm',
        'p4gm': '4mm',
        'p3': '3',
        'p3m1': '3m',
        'p31m': '3m',
        'p6': '6',
        'p6mm': '6mm'
    }

    # The order (number of symmetry operations) for each point group.
    point_group_orders = {
        '1': 1,
        '2': 2,
        'm': 2,
        '2mm': 4,
        '3': 3,
        '3m': 6,
        '4': 4,
        '4mm': 4,
        '6': 6,
        '6mm': 6
    }

    groups_with_pg_2 = []
    groups_with_pg_m = []

    # Iterate through all wallpaper groups to find those with a point group of order 2.
    for group, pg in wallpaper_to_point_group.items():
        if point_group_orders.get(pg) == 2:
            if pg == '2':
                groups_with_pg_2.append(group)
            elif pg == 'm':
                groups_with_pg_m.append(group)

    count_pg_2 = len(groups_with_pg_2)
    count_pg_m = len(groups_with_pg_m)
    total_count = count_pg_2 + count_pg_m

    print("The point groups of order 2 are '2' (rotation) and 'm' (reflection).")
    print("-" * 40)
    print(f"Number of wallpaper groups with point group '2': {count_pg_2}")
    print(f"They are: {', '.join(groups_with_pg_2)}\n")
    print(f"Number of wallpaper groups with point group 'm': {count_pg_m}")
    print(f"They are: {', '.join(groups_with_pg_m)}\n")
    print("To find the total, we sum the counts from each point group:")
    print(f"Equation: {count_pg_2} + {count_pg_m} = {total_count}")
    print("-" * 40)
    print(f"Total number of wallpaper groups with a point group of order 2 is: {total_count}")

solve_wallpaper_groups()