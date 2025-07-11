import collections

def solve_wallpaper_group_count():
    """
    This function finds and counts the number of wallpaper groups
    that have a point group of a specific order.
    """
    # Data for the 17 wallpaper groups. Each tuple contains:
    # (Wallpaper Group Name, Point Group Name, Order of the Point Group)
    wallpaper_groups_data = [
        ('p1', '1', 1),
        ('p2', '2', 2),
        ('pm', 'm', 2),
        ('pg', 'm', 2),
        ('cm', 'm', 2),
        ('p2mm', '2mm', 4),
        ('p2mg', '2mm', 4),
        ('p2gg', '2', 2),
        ('c2mm', '2mm', 4),
        ('p3', '3', 3),
        ('p3m1', '3m', 6),
        ('p31m', '3m', 6),
        ('p4', '4', 4),
        ('p4mm', '4mm', 8),
        ('p4gm', '4mm', 8),
        ('p6', '6', 6),
        ('p6mm', '6mm', 12),
    ]

    target_order = 2

    # Group the wallpaper groups by their point group name
    groups_by_point_group = collections.defaultdict(list)
    for group_name, point_group, order in wallpaper_groups_data:
        if order == target_order:
            groups_by_point_group[point_group].append(group_name)

    print(f"Finding wallpaper groups with a point group of order {target_order}:")
    print("-" * 30)

    counts = []
    total_count = 0

    # Sort keys for consistent output order
    for point_group in sorted(groups_by_point_group.keys()):
        groups = groups_by_point_group[point_group]
        count = len(groups)
        counts.append(str(count))
        total_count += count
        print(f"Point group '{point_group}' has {count} wallpaper group(s): {', '.join(groups)}")

    print("-" * 30)
    
    # Construct and print the final equation
    equation_str = " + ".join(counts)
    print(f"The total is the sum of these counts.")
    print(f"Final Equation: {equation_str} = {total_count}")


solve_wallpaper_group_count()
<<<5>>>