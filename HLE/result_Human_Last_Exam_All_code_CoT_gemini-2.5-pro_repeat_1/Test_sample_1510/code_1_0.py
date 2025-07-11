def count_wallpaper_groups():
    """
    This function calculates and prints the number of wallpaper groups
    that have a point group of order 2.
    """
    # A mapping of the 17 wallpaper groups to their respective point groups
    # using Hermann-Mauguin notation.
    wallpaper_groups_map = {
        'p1': '1',
        'p2': '2',
        'pm': 'm',
        'pg': 'm',
        'cm': 'm',
        'p2mm': '2mm',
        'p2mg': '2mm',
        'p2gg': '2mm',
        'c2mm': '2mm',
        'p4': '4',
        'p4mm': '4mm',
        'p4gm': '4mm',
        'p3': '3',
        'p3m1': '3m',
        'p31m': '3m',
        'p6': '6',
        'p6mm': '6mm'
    }

    # The point groups of order 2 are '2' and 'm'.
    point_groups_order_2 = ['2', 'm']

    # A dictionary to store the lists of wallpaper groups for each point group of order 2.
    result = {pg: [] for pg in point_groups_order_2}

    # Iterate through all wallpaper groups and categorize them.
    for group, point_group in wallpaper_groups_map.items():
        if point_group in point_groups_order_2:
            result[point_group].append(group)

    print("Finding the wallpaper groups with a point group of order 2...")
    print("-" * 60)

    total_count = 0
    equation_parts = []

    # Print the breakdown for each point group of order 2.
    for point_group, groups in result.items():
        count = len(groups)
        print(f"Point group '{point_group}' has {count} wallpaper group(s): {', '.join(groups)}")
        total_count += count
        equation_parts.append(str(count))

    print("-" * 60)
    # Print the final equation and the total count.
    final_equation = " + ".join(equation_parts)
    print(f"The total number is the sum of these counts.")
    print(f"Total = {final_equation} = {total_count}")

# Execute the function
count_wallpaper_groups()