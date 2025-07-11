def count_wallpaper_groups_by_point_group_order():
    """
    This function identifies and counts the wallpaper groups that have a point group of a specific order.
    """
    # A dictionary mapping each of the 17 wallpaper groups to its point group name and order.
    wallpaper_group_data = {
        "p1":   ("1", 1),
        "p2":   ("2", 2),
        "pm":   ("m", 2),
        "pg":   ("m", 2),
        "cm":   ("m", 2),
        "p2mm": ("2mm", 4),
        "p2mg": ("2mm", 4),
        "p2gg": ("2mm", 4),
        "c2mm": ("2mm", 4),
        "p4":   ("4", 4),
        "p4mm": ("4mm", 8),
        "p4gm": ("4mm", 8),
        "p3":   ("3", 3),
        "p3m1": ("3m", 6),
        "p31m": ("3m", 6),
        "p6":   ("6", 6),
        "p6mm": ("6mm", 12),
    }

    target_order = 2
    
    # Find groups with point group '2' (order 2)
    groups_with_pg_2 = [name for name, data in wallpaper_group_data.items() if data == ("2", target_order)]
    
    # Find groups with point group 'm' (order 2)
    groups_with_pg_m = [name for name, data in wallpaper_group_data.items() if data == ("m", target_order)]

    count_pg_2 = len(groups_with_pg_2)
    count_pg_m = len(groups_with_pg_m)
    total_count = count_pg_2 + count_pg_m
    
    print(f"There is {count_pg_2} wallpaper group with point group '2' (order 2): {', '.join(groups_with_pg_2)}")
    print(f"There are {count_pg_m} wallpaper groups with point group 'm' (order 2): {', '.join(groups_with_pg_m)}")
    print("\nTo find the total number of wallpaper groups with a point group of order 2, we sum the counts:")
    print(f"{count_pg_2} + {count_pg_m} = {total_count}")
    
count_wallpaper_groups_by_point_group_order()