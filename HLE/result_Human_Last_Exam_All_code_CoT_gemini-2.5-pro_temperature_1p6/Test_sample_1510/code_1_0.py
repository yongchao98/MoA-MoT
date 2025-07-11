def count_wallpaper_groups_by_point_group_order():
    """
    This function identifies and counts the wallpaper groups that have a point group of order 2.
    It prints a detailed breakdown of the groups and the final calculation.
    """
    # A dictionary mapping each of the 17 wallpaper groups to its
    # point group (crystal class) and the order of that point group.
    # The notation used is: Hermann-Mauguin for wallpaper groups, Schoenflies for point groups.
    # D1 is also commonly denoted as 'm'.
    wallpaper_groups_data = {
        'p1': ('C1', 1),
        'p2': ('C2', 2),
        'pm': ('D1', 2),
        'pg': ('D1', 2),
        'cm': ('D1', 2),
        'p2mm': ('D2', 4),
        'p2mg': ('D2', 4),
        'p2gg': ('D2', 4),
        'c2mm': ('D2', 4),
        'p4': ('C4', 4),
        'p4mm': ('D4', 8),
        'p4gm': ('D4', 8),
        'p3': ('C3', 3),
        'p3m1': ('D3', 6),
        'p31m': ('D3', 6),
        'p6': ('C6', 6),
        'p6mm': ('D6', 12)
    }

    target_order = 2

    # Group the wallpaper groups by their point group name if the order is 2
    groups_by_point_group = {}
    for group_name, (pg_name, pg_order) in wallpaper_groups_data.items():
        if pg_order == target_order:
            if pg_name not in groups_by_point_group:
                groups_by_point_group[pg_name] = []
            groups_by_point_group[pg_name].append(group_name)

    print(f"Finding wallpaper groups with a point group of order {target_order}:\n")

    counts = []
    # Print the detailed results
    for pg_name in sorted(groups_by_point_group.keys()):
        groups = groups_by_point_group[pg_name]
        count = len(groups)
        counts.append(count)
        print(f"Point Group {pg_name} (order {target_order}):")
        print(f"  - Count: {count}")
        print(f"  - Groups: {', '.join(sorted(groups))}")
        print("-" * 25)

    # Print the final calculation and total
    total_count = sum(counts)
    calculation_string = " + ".join(map(str, counts))
    
    print("\nFinal Calculation:")
    print(f"{calculation_string} = {total_count}")


if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()