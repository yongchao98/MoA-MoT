def count_wallpaper_groups():
    """
    This function identifies and counts the wallpaper groups
    that have a point group of order 2.
    """
    # The 2D crystallographic point groups of order 2 are C2 and D1 (also denoted as 'm').
    
    # List the wallpaper groups associated with each point group of order 2.
    groups_by_point_group = {
        "C2": ["p2"],
        "D1 (m)": ["pm", "pg", "cm"]
    }

    # Count the number of groups for each point group.
    count_c2 = len(groups_by_point_group["C2"])
    count_d1_m = len(groups_by_point_group["D1 (m)"])
    
    # Calculate the total.
    total_groups = count_c2 + count_d1_m

    print("The crystallographic point groups of order 2 are C2 and D1 (m).")
    print(f"\nNumber of wallpaper groups with point group C2: {count_c2}")
    print(f"Number of wallpaper groups with point group D1 (m): {count_d1_m}")

    # The final equation as requested.
    print(f"\nTotal number = {count_c2} (for C2) + {count_d1_m} (for D1) = {total_groups}")

count_wallpaper_groups()