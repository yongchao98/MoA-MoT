def solve_wallpaper_group_count():
    """
    This function identifies and counts the wallpaper groups
    that have a point group of order 2.
    """
    # A dictionary mapping each of the 17 wallpaper groups to the order of its point group.
    # The point groups of order 2 are C2 and D1 (m).
    wallpaper_group_data = {
        "p1": 1,
        "p2": 2,    # Point group C2
        "pm": 2,    # Point group D1 (m)
        "pg": 2,    # Point group D1 (m)
        "cm": 2,    # Point group D1 (m)
        "p2mm": 4,
        "p2mg": 4,
        "p2gg": 4,
        "c2mm": 4,
        "p4": 4,
        "p4mm": 8,
        "p4gm": 8,
        "p3": 3,
        "p3m1": 6,
        "p31m": 6,
        "p6": 6,
        "p6m": 12
    }

    # Filter the dictionary to find groups with a point group of order 2.
    groups_with_order_2 = []
    for group_name, order in wallpaper_group_data.items():
        if order == 2:
            groups_with_order_2.append(group_name)

    # The total count is the number of groups found.
    count = len(groups_with_order_2)

    print(f"The wallpaper groups with a point group of order 2 are: {', '.join(groups_with_order_2)}.")

    # Display the count as a sum, showing '1' for each group found.
    sum_equation = " + ".join(['1'] * count)
    print(f"The final calculation is: {sum_equation} = {count}")

solve_wallpaper_group_count()