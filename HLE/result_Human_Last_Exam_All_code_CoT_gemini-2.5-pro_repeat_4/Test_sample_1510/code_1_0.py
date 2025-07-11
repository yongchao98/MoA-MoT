def count_wallpaper_groups_by_order():
    """
    This script identifies and counts the wallpaper groups
    that have a point group of order 2.
    """

    # A mapping of the 17 wallpaper groups to their respective point groups.
    wallpaper_to_point_group = {
        "p1": "1",
        "p2": "2",
        "pm": "m",
        "pg": "m",
        "cm": "m",
        "p2mm": "2mm",
        "p2mg": "2mm",
        "p2gg": "2",
        "c2mm": "2mm",
        "p4": "4",
        "p4mm": "4mm",
        "p4gm": "4mm",
        "p3": "3",
        "p3m1": "3m",
        "p31m": "3m",
        "p6": "6",
        "p6mm": "6mm",
    }

    # The order for each of the 10 2D crystallographic point groups.
    # The point groups with order 2 are '2' and 'm'.
    point_group_orders = {
        '1': 1,
        '2': 2,    # Consists of {identity, 180-degree rotation}
        'm': 2,    # Consists of {identity, reflection}
        '2mm': 4,
        '3': 3,
        '3m': 6,
        '4': 4,
        '4mm': 8,
        '6': 6,
        '6mm': 12,
    }

    # Lists to hold the groups that match our criteria.
    groups_with_pg_2 = []
    groups_with_pg_m = []

    # Iterate through all wallpaper groups and categorize them.
    for group, point_group in wallpaper_to_point_group.items():
        if point_group_orders.get(point_group) == 2:
            if point_group == '2':
                groups_with_pg_2.append(group)
            elif point_group == 'm':
                groups_with_pg_m.append(group)

    print("The 2D crystallographic point groups of order 2 are '2' (rotational) and 'm' (reflectional).\n")

    print(f"There are {len(groups_with_pg_2)} wallpaper groups with point group '2':")
    print(", ".join(sorted(groups_with_pg_2)))
    print("-" * 30)

    print(f"There are {len(groups_with_pg_m)} wallpaper groups with point group 'm':")
    print(", ".join(sorted(groups_with_pg_m)))
    print("-" * 30)

    total = len(groups_with_pg_2) + len(groups_with_pg_m)
    print("The total number of wallpaper groups with a point group of order 2 is:")
    print(f"{len(groups_with_pg_2)} + {len(groups_with_pg_m)} = {total}")

count_wallpaper_groups_by_order()