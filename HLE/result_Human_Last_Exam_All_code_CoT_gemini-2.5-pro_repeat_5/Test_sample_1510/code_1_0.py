def count_wallpaper_groups_by_point_group_order():
    """
    This script identifies and counts the wallpaper groups whose point group has an order of 2.
    """
    # A dictionary mapping each of the 17 wallpaper groups to its point group.
    wallpaper_to_point_group = {
        'p1': 'C1',
        'p2': 'C2',
        'pm': 'D1',
        'pg': 'D1',
        'cm': 'D1',
        'p2mm': 'D2',
        'p2mg': 'D2',
        'p2gg': 'D2',
        'c2mm': 'D2',
        'p4': 'C4',
        'p4mm': 'D4',
        'p4gm': 'D4',
        'p3': 'C3',
        'p3m1': 'D3',
        'p31m': 'D3',
        'p6': 'C6',
        'p6mm': 'D6',
    }

    # The point groups of order 2 are C2 and D1.
    point_group_c2 = 'C2'
    point_group_d1 = 'D1'

    # Count how many wallpaper groups fall into each category.
    count_c2 = 0
    count_d1 = 0

    list_c2 = []
    list_d1 = []

    for group, point_group in wallpaper_to_point_group.items():
        if point_group == point_group_c2:
            count_c2 += 1
            list_c2.append(group)
        elif point_group == point_group_d1:
            count_d1 += 1
            list_d1.append(group)
    
    total_count = count_c2 + count_d1

    print(f"The point groups of order 2 are {point_group_c2} and {point_group_d1}.")
    print("-" * 30)
    print(f"Number of wallpaper groups with point group {point_group_c2}: {count_c2} (Group: {', '.join(list_c2)})")
    print(f"Number of wallpaper groups with point group {point_group_d1}: {count_d1} (Groups: {', '.join(list_d1)})")
    print("-" * 30)
    print("The total number of wallpaper groups with a point group of order 2 is calculated as:")
    print(f"{count_c2} + {count_d1} = {total_count}")

# Execute the function
count_wallpaper_groups_by_point_group_order()