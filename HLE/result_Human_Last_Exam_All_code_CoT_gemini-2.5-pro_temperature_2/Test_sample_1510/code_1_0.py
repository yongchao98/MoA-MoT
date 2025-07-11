def count_wallpaper_groups_by_point_group_order():
    """
    Identifies and counts the wallpaper groups with a specific point group order.
    The 17 wallpaper groups are pre-defined with their respective point group orders.
    """
    # Data for the 17 wallpaper groups.
    # Each entry contains the group name and the order of its associated point group.
    wallpaper_groups_data = {
        'p1': 1, 'p2': 2, 'pm': 2, 'pg': 2, 'cm': 2, 'pmm': 4,
        'pmg': 4, 'pgg': 4, 'cmm': 4, 'p4': 4, 'p4m': 8, 'p4g': 8,
        'p3': 3, 'p3m1': 6, 'p31m': 6, 'p6': 6, 'p6m': 12
    }
    
    target_order = 2
    
    # Find groups with the target point group order
    matching_groups = []
    for group_name, order in wallpaper_groups_data.items():
        if order == target_order:
            matching_groups.append(group_name)
    
    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in matching_groups:
        print(f"- {group}")
        
    count = len(matching_groups)
    
    # Constructing and printing the final equation as requested
    equation_parts = ['1'] * count
    equation_str = " + ".join(equation_parts)
    
    print(f"\nThere are {count} such groups.")
    print(f"The calculation is: {equation_str} = {count}")

count_wallpaper_groups_by_point_group_order()