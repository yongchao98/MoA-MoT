def count_wallpaper_groups():
    """
    This function counts the number of wallpaper groups with a point group of order 2.
    It stores all 17 wallpaper groups and their point group orders in a dictionary,
    then filters them to find the ones that match the condition.
    """
    # Data for the 17 wallpaper groups and their point group orders.
    # The key is the group name (IUC notation) and the value is the order of its point group.
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

    # Find the groups with a point group of order 2
    groups_of_order_2 = [name for name, order in wallpaper_group_data.items() if order == 2]

    # Print the results
    print(f"The wallpaper groups with a point group of order 2 are: {', '.join(groups_of_order_2)}")

    count = len(groups_of_order_2)
    
    # Create the equation string like "1 + 1 + ... = total"
    equation_parts = ['1'] * count
    equation_str = ' + '.join(equation_parts)

    print(f"There are a total of {count} such groups.")
    print(f"Equation: {equation_str} = {count}")


count_wallpaper_groups()
<<<4>>>