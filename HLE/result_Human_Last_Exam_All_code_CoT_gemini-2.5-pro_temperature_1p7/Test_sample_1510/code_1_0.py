def count_wallpaper_groups_by_order():
    """
    Counts the number of wallpaper groups with a point group of order 2.
    """
    # A dictionary mapping each wallpaper group to its point group and the order of that point group.
    wallpaper_groups_info = {
        'p1': ('C1', 1),
        'p2': ('C2', 2),
        'pm': ('D1', 2),
        'pg': ('D1', 2),
        'cm': ('D1', 2),
        'p2mm': ('D2', 4),
        'p2mg': ('D2', 4),
        'p2gg': ('D2', 4),
        'c2mm': ('D2', 4),
        'p3': ('C3', 3),
        'p3m1': ('D3', 6),
        'p31m': ('D3', 6),
        'p4': ('C4', 4),
        'p4mm': ('D4', 8),
        'p4gm': ('D4', 8),
        'p6': ('C6', 6),
        'p6mm': ('D6', 12)
    }

    groups_with_c2 = []
    groups_with_d1 = []

    # Iterate through the groups and categorize them
    for group, (point_group, order) in wallpaper_groups_info.items():
        if order == 2:
            if point_group == 'C2':
                groups_with_c2.append(group)
            elif point_group == 'D1':
                groups_with_d1.append(group)

    num_c2 = len(groups_with_c2)
    num_d1 = len(groups_with_d1)
    total_count = num_c2 + num_d1
    
    print(f"Wallpaper groups with point group C2 (order 2): {', '.join(groups_with_c2)}")
    print(f"Number of C2 groups: {num_c2}")
    print("-" * 30)
    print(f"Wallpaper groups with point group D1 (order 2): {', '.join(groups_with_d1)}")
    print(f"Number of D1 groups: {num_d1}")
    print("-" * 30)
    print(f"The total number of wallpaper groups with point group of order 2 is calculated by summing the counts from each type.")
    print(f"Calculation: {num_c2} (from C2) + {num_d1} (from D1) = {total_count}")

count_wallpaper_groups_by_order()