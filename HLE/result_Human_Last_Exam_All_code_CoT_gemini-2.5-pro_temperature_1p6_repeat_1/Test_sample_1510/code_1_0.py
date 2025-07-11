import collections

def count_wallpaper_groups_by_point_group_order():
    """
    Identifies and counts the wallpaper groups with a point group of a specific order.
    """
    # Data for the 17 wallpaper groups. Each entry contains the name,
    # the associated point group, and the order of that point group.
    wallpaper_groups_data = [
        ("p1", "C1", 1),
        ("p2", "C2", 2),
        ("pm", "D1", 2),
        ("pg", "D1", 2),
        ("cm", "D1", 2),
        ("p2mm", "D2", 4),
        ("p2mg", "D2", 4),
        ("p2gg", "D2", 4),
        ("c2mm", "D2", 4),
        ("p4", "C4", 4),
        ("p4mm", "D4", 8),
        ("p4gm", "D4", 8),
        ("p3", "C3", 3),
        ("p3m1", "D3", 6),
        ("p31m", "D3", 6),
        ("p6", "C6", 6),
        ("p6m", "D6", 12),
    ]

    target_order = 2
    
    # Use a dictionary to store lists of groups, keyed by their point group name
    found_groups = collections.defaultdict(list)

    for name, point_group, order in wallpaper_groups_data:
        if order == target_order:
            found_groups[point_group].append(name)
    
    total_count = 0
    equation_parts = []

    print(f"Finding wallpaper groups with a point group of order {target_order}:\n")

    # Sort the point groups for consistent output
    sorted_point_groups = sorted(found_groups.keys())

    for point_group in sorted_point_groups:
        groups = found_groups[point_group]
        count = len(groups)
        total_count += count
        equation_parts.append(str(count))
        print(f"- Point group {point_group} (order 2) corresponds to {count} wallpaper group(s): {', '.join(groups)}")

    # Print the final equation and the total count
    equation = " + ".join(equation_parts)
    print(f"\nThe total number is the sum of the counts for each point group type.")
    print(f"Final calculation: {equation} = {total_count}")


if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()