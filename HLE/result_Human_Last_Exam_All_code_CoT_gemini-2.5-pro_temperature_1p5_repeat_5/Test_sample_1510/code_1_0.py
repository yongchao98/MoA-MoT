def count_wallpaper_groups_by_point_group_order():
    """
    This function calculates and displays the number of wallpaper groups
    that have a point group of a specific order.
    """
    # A dictionary mapping each of the 17 wallpaper groups to the order
    # of its corresponding crystallographic point group.
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
    found_groups = []

    # Iterate through the dictionary to find groups with the target point group order
    for group, order in wallpaper_group_data.items():
        if order == target_order:
            found_groups.append(group)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in found_groups:
        print(f"- {group}")

    count = len(found_groups)
    
    # Create an equation string showing each group counted as '1'
    equation_parts = ['1'] * count
    equation_str = " + ".join(equation_parts)

    print(f"\nFinal count calculation: {equation_str} = {count}")


if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()