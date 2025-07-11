def count_wallpaper_groups():
    """
    This script identifies and counts the wallpaper groups whose point group has an order of 2.
    """
    # A dictionary mapping the 17 wallpaper groups (IUC notation) to the order of their point group.
    wallpaper_group_data = {
        # Group Name: Point Group Order
        'p1': 1,
        'p2': 2,
        'pm': 2,
        'pg': 2,
        'cm': 2,
        'pmm': 4,
        'pmg': 4,
        'pgg': 4,
        'cmm': 4,
        'p4': 4,
        'p4m': 8,
        'p4g': 8,
        'p3': 3,
        'p3m1': 6,
        'p31m': 6,
        'p6': 6,
        'p6m': 12
    }

    # Find the groups with a point group of order 2
    groups_with_order_2 = []
    for group, order in wallpaper_group_data.items():
        if order == 2:
            groups_with_order_2.append(group)

    # Print the names of the identified groups
    print(f"The wallpaper groups with a point group of order 2 are: {', '.join(groups_with_order_2)}")

    # Create and print the equation showing the count
    count = len(groups_with_order_2)
    equation_numbers = ['1'] * count
    equation_str = " + ".join(equation_numbers)
    
    print(f"The calculation for the total count is: {equation_str} = {count}")

# Execute the function
count_wallpaper_groups()