def count_wallpaper_groups_by_point_group_order():
    """
    This function determines how many of the 17 wallpaper groups
    have a crystallographic point group of order 2.
    """
    # Data for the 17 wallpaper groups, their point groups, and the order of those point groups.
    # Note: 'pg' and 'cm' have point group 'm', as glide reflections and mirror reflections
    # are equivalent in the point group.
    wallpaper_data = [
        {'name': 'p1', 'point_group': '1', 'order': 1},
        {'name': 'p2', 'point_group': '2', 'order': 2},
        {'name': 'pm', 'point_group': 'm', 'order': 2},
        {'name': 'pg', 'point_group': 'm', 'order': 2},
        {'name': 'cm', 'point_group': 'm', 'order': 2},
        {'name': 'p2mm', 'point_group': '2mm', 'order': 4},
        {'name': 'p2mg', 'point_group': '2mm', 'order': 4},
        {'name': 'p2gg', 'point_group': '2mm', 'order': 4},
        {'name': 'c2mm', 'point_group': '2mm', 'order': 4},
        {'name': 'p4', 'point_group': '4', 'order': 4},
        {'name': 'p4mm', 'point_group': '4mm', 'order': 8},
        {'name': 'p4gm', 'point_group': '4mm', 'order': 8},
        {'name': 'p3', 'point_group': '3', 'order': 3},
        {'name': 'p3m1', 'point_group': '3m', 'order': 6},
        {'name': 'p31m', 'point_group': '3m', 'order': 6},
        {'name': 'p6', 'point_group': '6', 'order': 6},
        {'name': 'p6mm', 'point_group': '6mm', 'order': 12},
    ]

    # Find the groups with point group order 2
    groups_with_order_2 = [
        group for group in wallpaper_data if group['order'] == 2
    ]

    # Separate them by their specific point group for clarity
    groups_with_pg_2 = [
        group['name'] for group in groups_with_order_2 if group['point_group'] == '2'
    ]
    groups_with_pg_m = [
        group['name'] for group in groups_with_order_2 if group['point_group'] == 'm'
    ]
    
    count_pg_2 = len(groups_with_pg_2)
    count_pg_m = len(groups_with_pg_m)
    total = count_pg_2 + count_pg_m
    
    print("The point groups of order 2 are '2' and 'm'.")
    print(f"\nWallpaper groups with point group '2': {', '.join(groups_with_pg_2)}")
    print(f"Count: {count_pg_2}")
    
    print(f"\nWallpaper groups with point group 'm': {', '.join(groups_with_pg_m)}")
    print(f"Count: {count_pg_m}")

    print("\nThe total number of wallpaper groups with a point group of order 2 is the sum:")
    # Print the final equation with the numbers
    print(f"{count_pg_2} + {count_pg_m} = {total}")

if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()