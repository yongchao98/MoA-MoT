def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of wallpaper groups that have a point group of order 2.
    """

    # Step 1: Define the 17 wallpaper groups and their corresponding point groups.
    # Note: D1 is also known as Cs or C1v. D2 is C2v, D3 is C3v, etc.
    wallpaper_to_point_group = {
        'p1': 'C1',
        'p2': 'C2',
        'pm': 'D1',
        'pg': 'D1',
        'cm': 'D1',
        'p2mm': 'D2',
        'p2mg': 'D2',
        'p2gg': 'D2',
        'c2mm': 'D2',
        'p3': 'C3',
        'p3m1': 'D3',
        'p31m': 'D3',
        'p4': 'C4',
        'p4mm': 'D4',
        'p4gm': 'D4',
        'p6': 'C6',
        'p6mm': 'D6'
    }

    # Step 2: Define the orders of the relevant point groups.
    # The order of a group is the number of its elements (symmetry operations).
    point_group_orders = {
        'C1': 1, 'C2': 2, 'D1': 2, 'D2': 4,
        'C3': 3, 'D3': 6, 'C4': 4, 'D4': 8,
        'C6': 6, 'D6': 12
    }

    # Step 3: Find the wallpaper groups whose point group has an order of 2.
    # We will categorize them by their specific point group.
    
    # Point group C2 (2-fold rotation)
    groups_c2 = [group for group, pg in wallpaper_to_point_group.items()
                 if pg == 'C2' and point_group_orders[pg] == 2]
    
    # Point group D1 (one reflection axis, also called Cs)
    groups_d1 = [group for group, pg in wallpaper_to_point_group.items()
                 if pg == 'D1' and point_group_orders[pg] == 2]

    count_c2 = len(groups_c2)
    count_d1 = len(groups_d1)
    total_count = count_c2 + count_d1

    print("The point groups of order 2 are C2 (2-fold rotation) and D1 (reflection).\n")

    print(f"Number of wallpaper groups with point group C2: {count_c2}")
    print(f"Groups: {', '.join(groups_c2)}\n")

    print(f"Number of wallpaper groups with point group D1: {count_d1}")
    print(f"Groups: {', '.join(groups_d1)}\n")
    
    print("Total number of wallpaper groups with a point group of order 2:")
    # Final equation with each number printed explicitly
    print(f"{count_c2} + {count_d1} = {total_count}")


if __name__ == '__main__':
    count_wallpaper_groups_by_point_group_order()