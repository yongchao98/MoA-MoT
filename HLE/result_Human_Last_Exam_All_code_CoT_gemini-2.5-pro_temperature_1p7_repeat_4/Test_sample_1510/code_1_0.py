def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups with a point group of a specific order.
    """
    # Data mapping the 17 wallpaper groups to their associated point group and the order of that point group.
    # The point groups for 2D crystallography are Cn and Dn.
    wallpaper_group_data = {
        'p1':   {'point_group': 'C1', 'order': 1},
        'p2':   {'point_group': 'C2', 'order': 2},
        'pm':   {'point_group': 'D1', 'order': 2},
        'pg':   {'point_group': 'D1', 'order': 2},
        'cm':   {'point_group': 'D1', 'order': 2},
        'pmm':  {'point_group': 'D2', 'order': 4},
        'pmg':  {'point_group': 'D2', 'order': 4},
        'pgg':  {'point_group': 'D2', 'order': 4},
        'cmm':  {'point_group': 'D2', 'order': 4},
        'p4':   {'point_group': 'C4', 'order': 4},
        'p4m':  {'point_group': 'D4', 'order': 8},
        'p4g':  {'point_group': 'D4', 'order': 8},
        'p3':   {'point_group': 'C3', 'order': 3},
        'p3m1': {'point_group': 'D3', 'order': 6},
        'p31m': {'point_group': 'D3', 'order': 6},
        'p6':   {'point_group': 'C6', 'order': 6},
        'p6m':  {'point_group': 'D6', 'order': 12}
    }

    target_order = 2

    # Find the point groups that have the target order
    point_groups_with_target_order = {
        data['point_group'] for data in wallpaper_group_data.values() if data['order'] == target_order
    }

    print(f"The point groups of order {target_order} are: {', '.join(sorted(list(point_groups_with_target_order)))}\n")

    # Group wallpaper groups by their point group for detailed output
    groups_by_pg = {}
    for name, data in wallpaper_group_data.items():
        if data['order'] == target_order:
            pg = data['point_group']
            if pg not in groups_by_pg:
                groups_by_pg[pg] = []
            groups_by_pg[pg].append(name)
            
    # Print the breakdown and prepare for the final equation
    counts = []
    print("Finding the number of groups for each point group type:")
    for pg in sorted(groups_by_pg.keys()):
        count = len(groups_by_pg[pg])
        counts.append(count)
        group_list = ", ".join(groups_by_pg[pg])
        print(f"- Point group {pg} has {count} wallpaper group(s): {group_list}")

    # Calculate and print the final equation
    total_count = sum(counts)
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    
    print("\nThe final equation is:")
    print(f"{equation_str} = {total_count}")

if __name__ == '__main__':
    count_wallpaper_groups_by_point_group_order()