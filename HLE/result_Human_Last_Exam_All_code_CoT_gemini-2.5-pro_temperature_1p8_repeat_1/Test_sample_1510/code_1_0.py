def solve_wallpaper_group_question():
    """
    This function determines and counts the number of wallpaper groups
    that have a point group of order 2.
    """
    # Step 1: Define all 17 wallpaper groups and their associated point groups.
    # The notation for wallpaper groups is Hermann-Mauguin (e.g., 'p1', 'p2mm').
    # The notation for point groups is Schoenflies (e.g., 'C1', 'D2').
    wallpaper_groups_data = {
        'p1': 'C1',
        'p2': 'C2',
        'pm': 'D1',
        'pg': 'D1',
        'cm': 'D1',
        'pmm': 'D2',
        'pmg': 'D2',
        'pgg': 'D2',
        'cmm': 'D2',
        'p3': 'C3',
        'p3m1': 'D3',
        'p31m': 'D3',
        'p4': 'C4',
        'p4m': 'D4',
        'p4g': 'D4',
        'p6': 'C6',
        'p6m': 'D6',
    }

    # Step 2: The point groups of order 2 are C2 and D1.
    # C2 represents a 2-fold rotation.
    # D1 represents a single reflection.

    # Step 3: Filter the wallpaper groups.
    groups_with_c2 = []
    groups_with_d1 = []

    for group_name, point_group in wallpaper_groups_data.items():
        if point_group == 'C2':
            groups_with_c2.append(group_name)
        elif point_group == 'D1':
            groups_with_d1.append(group_name)

    # Step 4: Print the results and the final calculation.
    print("The point groups of order 2 are C2 and D1.")
    print("-" * 40)

    print(f"Number of wallpaper groups with point group C2: {len(groups_with_c2)}")
    print(f"Group(s): {', '.join(groups_with_c2)}")
    print("-" * 40)

    print(f"Number of wallpaper groups with point group D1: {len(groups_with_d1)}")
    print(f"Group(s): {', '.join(groups_with_d1)}")
    print("-" * 40)

    total_count = len(groups_with_c2) + len(groups_with_d1)
    print("The total number of wallpaper groups with point group of order 2 is the sum:")
    # The final print statement outputs each number in the equation as requested.
    print(f"{len(groups_with_c2)} + {len(groups_with_d1)} = {total_count}")

# Execute the function to solve the problem.
solve_wallpaper_group_question()