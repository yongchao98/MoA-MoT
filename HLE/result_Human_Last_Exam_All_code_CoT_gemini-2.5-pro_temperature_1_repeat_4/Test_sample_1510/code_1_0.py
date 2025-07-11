def find_wallpaper_groups_by_order():
    """
    This function identifies and counts the wallpaper groups with a point group of a specific order.
    """
    # A dictionary mapping each of the 17 wallpaper groups to the order of its point group.
    wallpaper_group_data = {
        'p1': 1,
        'p2': 2,
        'pm': 2,
        'pg': 2,
        'cm': 2,
        'p2mm': 4,
        'p2mg': 4,
        'p2gg': 4,
        'c2mm': 4,
        'p4': 4,
        'p4mm': 8,
        'p4gm': 8,
        'p3': 3,
        'p3m1': 6,
        'p31m': 6,
        'p6': 6,
        'p6mm': 12
    }

    target_order = 2
    
    # Find the groups with a point group of the target order.
    found_groups = [name for name, order in wallpaper_group_data.items() if order == target_order]
    
    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group_name in found_groups:
        print(f"- {group_name}")
    
    # Build and print the final equation.
    count = len(found_groups)
    if count > 0:
        equation_parts = ['1'] * count
        equation_str = " + ".join(equation_parts)
        print(f"\nThe total number is: {equation_str} = {count}")
    else:
        print(f"\nThere are 0 groups with a point group of order {target_order}.")

# Execute the function to find the answer.
find_wallpaper_groups_by_order()