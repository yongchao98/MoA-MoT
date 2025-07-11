def count_wallpaper_groups_by_point_group_order():
    """
    This function identifies and counts the wallpaper groups that have a point group of a specific order.
    It then prints the names of these groups and the final count as a sum.
    """
    # The 17 wallpaper groups and the order of their corresponding point group.
    # Data is in the format: ('group_name', point_group_order)
    wallpaper_groups_data = {
        'p1': 1,
        'p2': 2,
        'pm': 2,
        'pg': 2,
        'cm': 2,
        'p2mm': 4,
        'p2mg': 4,
        'p2gg': 4,
        'c2mm': 4,
        'p3': 3,
        'p3m1': 6,
        'p31m': 6,
        'p4': 4,
        'p4mm': 8,
        'p4gm': 8,
        'p6': 6,
        'p6mm': 12
    }

    target_order = 2
    
    # Find the wallpaper groups with a point group of order 2
    found_groups = []
    for group_name, order in wallpaper_groups_data.items():
        if order == target_order:
            found_groups.append(group_name)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in found_groups:
        print(f"- {group}")
    
    count = len(found_groups)
    
    # Create the summation string as requested
    sum_equation = " + ".join(['1'] * count)

    print(f"\nIn total, there are {count} such groups.")
    print(f"This can be represented as the sum: {sum_equation} = {count}")

# Execute the function
count_wallpaper_groups_by_point_group_order()