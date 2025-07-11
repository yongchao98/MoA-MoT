import collections

def count_wallpaper_groups_by_point_group_order():
    """
    Identifies and counts the wallpaper groups that have a point group of a specific order.
    """
    # Step 1: Define the mapping from each wallpaper group to its point group.
    # The keys are the 17 wallpaper groups in IUC notation.
    # The values are their corresponding point groups in Hermann-Mauguin notation.
    wallpaper_to_point_group = {
        'p1': '1', 'p2': '2', 'pm': 'm', 'pg': 'm', 'cm': 'm',
        'pmm': '2mm', 'pmg': '2mm', 'pgg': '2mm', 'cmm': '2mm',
        'p4': '4', 'p4m': '4mm', 'p4g': '4mm',
        'p3': '3', 'p3m1': '3m', 'p31m': '3m',
        'p6': '6', 'p6m': '6mm'
    }

    # Step 2: Define the order of each point group.
    # The keys are the point groups, and values are their orders.
    point_group_orders = {
        '1': 1, '2': 2, 'm': 2, '2mm': 4, '4': 4, '4mm': 8,
        '3': 3, '3m': 6, '6': 6, '6mm': 12
    }
    
    target_order = 2
    found_groups = []

    # Step 3: Iterate through wallpaper groups to find those with a point group of the target order.
    for group, point_group in wallpaper_to_point_group.items():
        if point_group_orders.get(point_group) == target_order:
            found_groups.append(group)
    
    # Sort for consistent, readable output
    found_groups.sort()

    print(f"The wallpaper groups with point group of order {target_order} are: {', '.join(found_groups)}.")

    # Step 4: Display the final count as an equation.
    count = len(found_groups)
    if count > 0:
        equation_parts = ["1"] * count
        equation_str = " + ".join(equation_parts)
        print(f"The total count is {equation_str} = {count}.")
    else:
        print("There are 0 such groups.")

count_wallpaper_groups_by_point_group_order()