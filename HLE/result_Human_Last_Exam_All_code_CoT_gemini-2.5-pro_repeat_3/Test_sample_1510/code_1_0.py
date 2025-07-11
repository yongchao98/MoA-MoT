def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups that have a point group of a specific order.
    In this case, it finds the groups with a point group of order 2.
    """
    # A dictionary mapping each of the 17 wallpaper groups to the order of its point group.
    wallpaper_group_data = {
        # Oblique lattice
        'p1': 1,  # Point group C1
        'p2': 2,  # Point group C2
        # Rectangular lattice (primitive)
        'pm': 2,  # Point group D1 (cm)
        'pg': 2,  # Point group D1 (cm)
        'p2mm': 4, # Point group D2 (c2v)
        'p2mg': 4, # Point group D2 (c2v)
        'p2gg': 4, # Point group D2 (c2v)
        # Rectangular lattice (centered)
        'cm': 2,  # Point group D1 (cm)
        'c2mm': 4, # Point group D2 (c2v)
        # Square lattice
        'p4': 4,  # Point group C4
        'p4m': 8,  # Point group D4 (c4v)
        'p4g': 8,  # Point group D4 (c4v)
        # Hexagonal lattice
        'p3': 3,  # Point group C3
        'p3m1': 6, # Point group D3 (c3v)
        'p31m': 6, # Point group D3 (c3v)
        'p6': 6,  # Point group C6
        'p6m': 12, # Point group D6 (c6v)
    }

    target_order = 2
    found_groups = []
    
    for group_name, order in wallpaper_group_data.items():
        if order == target_order:
            found_groups.append(group_name)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in found_groups:
        print(group)

    count = len(found_groups)
    
    # Create the equation string
    # e.g., "1 + 1 + 1 + 1 = 4"
    equation_parts = ["1"] * count
    equation_str = " + ".join(equation_parts)
    
    print(f"\nThere are a total of {count} such groups.")
    if count > 0:
        print(f"The calculation is: {equation_str} = {count}")

count_wallpaper_groups_by_point_group_order()