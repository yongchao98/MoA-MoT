def solve_wallpaper_group_question():
    """
    This script finds and counts the number of wallpaper groups that have a point group of order 2.
    """
    # A dictionary mapping each of the 17 wallpaper groups to its point group and the order of that point group.
    wallpaper_groups_info = {
        # Group Name: (Point Group Symbol, Point Group Order)
        "p1":   ("C1", 1),
        "p2":   ("C2", 2), # Has a point group of order 2
        "pm":   ("Cm", 2), # Has a point group of order 2
        "pg":   ("Cm", 2), # Has a point group of order 2
        "cm":   ("Cm", 2), # Has a point group of order 2
        "p2mm": ("D2", 4),
        "p2mg": ("D2", 4),
        "p2gg": ("D2", 4),
        "c2mm": ("D2", 4),
        "p4":   ("C4", 4),
        "p4mm": ("D4", 8),
        "p4gm": ("D4", 8),
        "p3":   ("C3", 3),
        "p3m1": ("D3", 6),
        "p31m": ("D3", 6),
        "p6":   ("C6", 6),
        "p6mm": ("D6", 12),
    }

    # Find the groups with point group of order 2
    groups_with_c2 = []
    groups_with_cm = []

    for group_name, (point_group, order) in wallpaper_groups_info.items():
        if order == 2:
            if point_group == "C2":
                groups_with_c2.append(group_name)
            elif point_group == "Cm":
                groups_with_cm.append(group_name)

    count_c2 = len(groups_with_c2)
    count_cm = len(groups_with_cm)
    total_count = count_c2 + count_cm

    print("The point groups of order 2 are C2 (2-fold rotation) and Cm (reflection).")
    print("\nWallpaper groups with point group C2:")
    for group in groups_with_c2:
        print(f"- {group}")
    print(f"Count: {count_c2}")

    print("\nWallpaper groups with point group Cm:")
    for group in groups_with_cm:
        print(f"- {group}")
    print(f"Count: {count_cm}")

    print("\nTo find the total number of wallpaper groups with a point group of order 2, we sum the counts:")
    print(f"{count_c2} (for C2) + {count_cm} (for Cm) = {total_count}")

solve_wallpaper_group_question()