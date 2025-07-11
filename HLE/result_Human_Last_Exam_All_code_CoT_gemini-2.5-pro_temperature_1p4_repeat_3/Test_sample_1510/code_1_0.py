import collections

def count_wallpaper_groups_by_point_group_order():
    """
    This function calculates and prints the number of wallpaper groups
    with a point group of a specific order.
    """
    # The 17 wallpaper groups and the order of their respective point groups.
    # The point group (or crystal class) describes the symmetries (rotations, reflections)
    # that leave at least one point fixed.
    wallpaper_group_data = {
        'p1': 1,
        'p2': 2,
        'pm': 2,
        'pg': 2,
        'cm': 2,
        'pmm': 4,
        'pmg': 4,
        'pgg': 4,
        'cmm': 4,
        'p4': 4,
        'p4m': 8,
        'p4g': 8,
        'p3': 3,
        'p3m1': 6,
        'p31m': 6,
        'p6': 6,
        'p6m': 12
    }

    target_order = 2
    
    # Find the groups with a point group of the target order
    matching_groups = []
    for group_name, order in wallpaper_group_data.items():
        if order == target_order:
            matching_groups.append(group_name)

    # Sort the list for consistent output
    matching_groups.sort()

    count = len(matching_groups)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    print(", ".join(matching_groups))
    
    # Create the equation string as requested
    equation_str = " + ".join(['1'] * count)
    
    print("\nEach of these groups contributes 1 to the total count.")
    print(f"The final calculation is: {equation_str} = {count}")
    
count_wallpaper_groups_by_point_group_order()