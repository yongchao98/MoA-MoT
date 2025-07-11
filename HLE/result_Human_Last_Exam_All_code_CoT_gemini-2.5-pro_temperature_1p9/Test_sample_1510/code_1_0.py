def find_wallpaper_groups_by_point_group_order():
    """
    This function finds and counts the wallpaper groups with a point group of a specific order.
    """
    # Data for the 17 wallpaper groups, including their point group and the order of the point group.
    # The point group D1 is also known as Cs or C1h.
    wallpaper_groups_data = [
        {'name': 'p1',   'point_group': 'C1', 'order': 1},
        {'name': 'p2',   'point_group': 'C2', 'order': 2},
        {'name': 'pm',   'point_group': 'D1', 'order': 2},
        {'name': 'pg',   'point_group': 'D1', 'order': 2},
        {'name': 'cm',   'point_group': 'D1', 'order': 2},
        {'name': 'p2mm', 'point_group': 'D2', 'order': 4},
        {'name': 'p2mg', 'point_group': 'D2', 'order': 4},
        {'name': 'p2gg', 'point_group': 'D2', 'order': 4},
        {'name': 'c2mm', 'point_group': 'D2', 'order': 4},
        {'name': 'p3',   'point_group': 'C3', 'order': 3},
        {'name': 'p3m1', 'point_group': 'D3', 'order': 6},
        {'name': 'p31m', 'point_group': 'D3', 'order': 6},
        {'name': 'p4',   'point_group': 'C4', 'order': 4},
        {'name': 'p4mm', 'point_group': 'D4', 'order': 8},
        {'name': 'p4gm', 'point_group': 'D4', 'order': 8},
        {'name': 'p6',   'point_group': 'C6', 'order': 6},
        {'name': 'p6mm', 'point_group': 'D6', 'order': 12},
    ]

    target_order = 2
    
    print(f"The point groups of order {target_order} are C2 and D1 (Cs).")
    
    # Separate the groups by point group type for a clear explanation.
    groups_C2 = [g['name'] for g in wallpaper_groups_data if g['point_group'] == 'C2' and g['order'] == target_order]
    groups_D1 = [g['name'] for g in wallpaper_groups_data if g['point_group'] == 'D1' and g['order'] == target_order]
    
    count_C2 = len(groups_C2)
    count_D1 = len(groups_D1)
    total_count = count_C2 + count_D1
    
    print(f"\nNumber of wallpaper groups with point group C2: {count_C2}")
    print(f"Group(s): {', '.join(groups_C2)}")
    
    print(f"\nNumber of wallpaper groups with point group D1: {count_D1}")
    print(f"Group(s): {', '.join(groups_D1)}")
    
    print("\nThe total number of wallpaper groups with a point group of order 2 is:")
    # The final equation with each number printed explicitly.
    print(f"{count_C2} + {count_D1} = {total_count}")

if __name__ == "__main__":
    find_wallpaper_groups_by_point_group_order()