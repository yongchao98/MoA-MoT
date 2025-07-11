def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups with a point group of a specific order.
    """
    # Data for the 17 wallpaper groups. Each tuple contains:
    # (Group Name, Point Group Name, Order of the Point Group)
    wallpaper_groups_data = [
        ('p1', '1', 1),
        ('p2', '2', 2),
        ('pm', 'm', 2),
        ('pg', 'm', 2),
        ('cm', 'm', 2),
        ('p2mm', '2mm', 4),
        ('p2mg', '2mm', 4),
        ('p2gg', '2mm', 4),
        ('c2mm', '2mm', 4),
        ('p4', '4', 4),
        ('p4mm', '4mm', 8),
        ('p4gm', '4mm', 8),
        ('p3', '3', 3),
        ('p3m1', '3m', 6),
        ('p31m', '3m', 6),
        ('p6', '6', 6),
        ('p6m', '6mm', 12)
    ]

    target_order = 2
    
    # Find all groups with the target point group order
    matching_groups = []
    for group_name, point_group, order in wallpaper_groups_data:
        if order == target_order:
            matching_groups.append(group_name)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    print(", ".join(matching_groups))
    
    # Create the equation string
    equation_parts = ["1"] * len(matching_groups)
    equation_str = " + ".join(equation_parts)
    
    total_count = len(matching_groups)
    
    print(f"\nThere are {equation_str} = {total_count} such groups.")

count_wallpaper_groups_by_point_group_order()