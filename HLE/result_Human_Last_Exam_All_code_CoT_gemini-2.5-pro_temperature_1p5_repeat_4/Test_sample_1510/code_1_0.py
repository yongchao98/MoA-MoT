def count_wallpaper_groups_by_point_group_order():
    """
    Counts the number of wallpaper groups with a point group of a specific order.
    """
    # A dictionary mapping each wallpaper group to the order of its point group.
    wallpaper_group_point_group_orders = {
        # Oblique lattice
        'p1': 1, 'p2': 2,
        # Rectangular lattice
        'pm': 2, 'pg': 2, 'cm': 2,
        'p2mm': 4, 'p2mg': 4, 'p2gg': 4, 'c2mm': 4,
        # Square lattice
        'p4': 4, 'p4mm': 8, 'p4gm': 8,
        # Hexagonal lattice
        'p3': 3, 'p3m1': 6, 'p31m': 6,
        'p6': 6, 'p6mm': 12
    }

    target_order = 2
    found_groups = []

    # Find all groups with a point group of the target order
    for group, order in wallpaper_group_point_group_orders.items():
        if order == target_order:
            found_groups.append(group)
    
    count = len(found_groups)
    
    print(f"The wallpaper groups with a point group of order {target_order} are: {', '.join(found_groups)}.")
    
    # Create the equation string as requested
    if count > 0:
        equation_parts = ["1" for _ in found_groups]
        equation_str = " + ".join(equation_parts)
        print(f"The final count is derived from the equation: {equation_str} = {count}")
    else:
        print("No wallpaper groups found with the specified point group order.")

    print(f"Total number of wallpaper groups with point group of order {target_order} is: {count}")

count_wallpaper_groups_by_point_group_order()