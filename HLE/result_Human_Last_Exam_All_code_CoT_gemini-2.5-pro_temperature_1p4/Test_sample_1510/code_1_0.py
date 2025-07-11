def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups that have a point group
    of a specific order and prints the details.
    """
    # A dictionary mapping each wallpaper group to its point group and the order of that point group.
    # Point group D1 is often just denoted by 'm' for mirror.
    wallpaper_groups_info = {
        "p1":   {"point_group": "C1", "order": 1},
        "p2":   {"point_group": "C2", "order": 2},
        "pm":   {"point_group": "D1", "order": 2},
        "pg":   {"point_group": "D1", "order": 2},
        "cm":   {"point_group": "D1", "order": 2},
        "p2mm": {"point_group": "D2", "order": 4},
        "p2mg": {"point_group": "D2", "order": 4},
        "p2gg": {"point_group": "D2", "order": 4},
        "c2mm": {"point_group": "D2", "order": 4},
        "p3":   {"point_group": "C3", "order": 3},
        "p3m1": {"point_group": "D3", "order": 6},
        "p31m": {"point_group": "D3", "order": 6},
        "p4":   {"point_group": "C4", "order": 4},
        "p4mm": {"point_group": "D4", "order": 8},
        "p4gm": {"point_group": "D4", "order": 8},
        "p6":   {"point_group": "C6", "order": 6},
        "p6mm": {"point_group": "D6", "order": 12},
    }

    target_order = 2
    
    # Group the wallpaper groups by their point group if their order matches the target.
    groups_by_point_group = {}
    for group_name, data in wallpaper_groups_info.items():
        if data["order"] == target_order:
            point_group = data["point_group"]
            if point_group not in groups_by_point_group:
                groups_by_point_group[point_group] = []
            groups_by_point_group[point_group].append(group_name)
    
    print(f"Finding wallpaper groups with a point group of order {target_order}:")
    
    counts = []
    # Print the breakdown of groups found for each relevant point group.
    for point_group, group_list in sorted(groups_by_point_group.items()):
        count = len(group_list)
        counts.append(count)
        print(f" - Point group {point_group} has {count} wallpaper group(s): {', '.join(group_list)}")
        
    total_count = sum(counts)
    
    print(f"\nThe total number is the sum of these counts.")
    print("Final equation:")
    
    # Construct and print the final equation.
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {total_count}")

# Execute the function
count_wallpaper_groups_by_point_group_order()