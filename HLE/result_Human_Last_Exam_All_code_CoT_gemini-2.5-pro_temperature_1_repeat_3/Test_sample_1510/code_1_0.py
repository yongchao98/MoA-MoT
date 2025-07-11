def count_wallpaper_groups_by_point_group_order():
    """
    Finds and counts the wallpaper groups with a point group of a specific order.
    """
    # Data for the 17 wallpaper groups.
    # The dictionary maps the group name to the order of its point group.
    wallpaper_groups_data = {
        # Oblique lattice
        'p1': 1, 'p2': 2,
        # Rectangular lattice
        'pm': 2, 'pg': 2, 'cm': 2, 'p2mm': 4, 'p2mg': 4, 'p2gg': 4, 'c2mm': 4,
        # Hexagonal lattice
        'p3': 3, 'p3m1': 6, 'p31m': 6, 'p6': 6, 'p6mm': 12,
        # Square lattice
        'p4': 4, 'p4mm': 8, 'p4gm': 8
    }

    target_order = 2
    found_groups = []

    # Find groups with the target point group order
    for group_name, order in wallpaper_groups_data.items():
        if order == target_order:
            found_groups.append(group_name)

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in found_groups:
        print(f"- {group}")

    # Calculate and display the total count as an equation
    count = len(found_groups)
    if count > 0:
        # Creates a string like "1 + 1 + 1 + 1"
        equation_parts = ['1'] * count
        equation_str = " + ".join(equation_parts)
        print(f"\nThe total count is calculated as: {equation_str} = {count}")
    else:
        print("\nNo groups found with the specified point group order.")

    print(f"\nThus, there are {count} wallpaper groups with a point group of order 2.")


if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()