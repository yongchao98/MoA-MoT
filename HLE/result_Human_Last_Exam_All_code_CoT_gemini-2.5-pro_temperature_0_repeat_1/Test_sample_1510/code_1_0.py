def count_wallpaper_groups_by_point_group_order():
    """
    This function identifies and counts the wallpaper groups
    whose point group has a specific order.
    """
    # Step 1: Define the 17 wallpaper groups and their point group orders.
    # The data is stored as a dictionary where the key is the wallpaper group name
    # and the value is the order of its associated point group.
    wallpaper_groups_data = {
        'p1': 1, 'p2': 2, 'pm': 2, 'pg': 2, 'cm': 2, 'p2mm': 4,
        'p2mg': 4, 'p2gg': 4, 'c2mm': 4, 'p4': 4, 'p4mm': 8,
        'p4gm': 8, 'p3': 3, 'p3m1': 6, 'p31m': 6, 'p6': 6, 'p6mm': 12
    }

    target_order = 2
    
    # Step 2: Filter the groups to find those with a point group of order 2.
    found_groups = []
    for group, order in wallpaper_groups_data.items():
        if order == target_order:
            found_groups.append(group)
            
    # Step 3: Count and display the results.
    count = len(found_groups)
    
    print(f"The wallpaper groups with a point group of order {target_order} are: {', '.join(found_groups)}.")
    
    # Create the equation string as requested
    equation_parts = ['1'] * count
    equation_str = " + ".join(equation_parts)
    
    print(f"The total number is {equation_str} = {count}.")

# Execute the function
count_wallpaper_groups_by_point_group_order()