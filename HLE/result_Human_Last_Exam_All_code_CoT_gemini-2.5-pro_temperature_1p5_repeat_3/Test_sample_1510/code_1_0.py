def count_wallpaper_groups_by_point_group_order():
    """
    This script identifies and counts the wallpaper groups that have a point group of a specific order.
    """
    # The 17 wallpaper groups and the order of their corresponding point groups.
    # Data is represented as a dictionary where the key is the group's IUC notation
    # and the value is the order of its point group.
    wallpaper_group_data = {
        'p1': 1,   # Point group 1 (C1)
        'p2': 2,   # Point group 2 (C2)
        'pm': 2,   # Point group m (Cs)
        'pg': 2,   # Point group m (Cs)
        'cm': 2,   # Point group m (Cs)
        'pmm': 4,  # Point group 2mm (D2)
        'pmg': 4,  # Point group 2mm (D2)
        'pgg': 4,  # Point group 2mm (D2)
        'cmm': 4,  # Point group 2mm (D2)
        'p4': 4,   # Point group 4 (C4)
        'p4m': 8,  # Point group 4mm (D4)
        'p4g': 8,  # Point group 4mm (D4)
        'p3': 3,   # Point group 3 (C3)
        'p3m1': 6, # Point group 3m (D3)
        'p31m': 6, # Point group 3m (D3)
        'p6': 6,   # Point group 6 (C6)
        'p6m': 12  # Point group 6mm (D6)
    }

    target_order = 2
    found_groups = []

    # Find all groups with a point group of the target order
    for group_name, order in wallpaper_group_data.items():
        if order == target_order:
            found_groups.append(group_name)
    
    # Sort the list alphabetically for consistent output
    found_groups.sort()

    print(f"The wallpaper groups with a point group of order {target_order} are:")
    for group in found_groups:
        print(f"- {group}")

    total_count = len(found_groups)
    
    if total_count > 0:
        # Create the equation string, e.g., "1 + 1 + 1 + 1"
        equation = " + ".join(["1"] * total_count)
        print(f"\nThe total number is calculated as: {equation} = {total_count}")
    else:
        print("\nNo such groups found.")

count_wallpaper_groups_by_point_group_order()