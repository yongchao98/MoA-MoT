def count_wallpaper_groups_by_point_group_order():
    """
    This function counts the number of 2D wallpaper groups whose point group has an order of 2.
    """
    # A mapping of the 17 wallpaper groups (full international notation) to their point groups.
    wallpaper_group_data = {
        'p1': '1',
        'p2': '2',      # Point group of order 2
        'pm': 'm',      # Point group of order 2
        'pg': 'm',      # Point group of order 2
        'cm': 'm',      # Point group of order 2
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
        'p6mm': '6m'
    }

    # Find groups with point group '2' or 'm'
    groups_pg_2 = [name for name, pg in wallpaper_group_data.items() if pg == '2']
    groups_pg_m = [name for name, pg in wallpaper_group_data.items() if pg == 'm']

    count_pg_2 = len(groups_pg_2)
    count_pg_m = len(groups_pg_m)
    total_count = count_pg_2 + count_pg_m

    print("The point groups of order 2 are '2' (180-degree rotation) and 'm' (reflection).")
    print("\nCounting the wallpaper groups for each:")
    
    print(f"\nNumber of groups with point group '2': {count_pg_2}")
    print(f"  (Group: {', '.join(groups_pg_2)})")
    
    print(f"Number of groups with point group 'm': {count_pg_m}")
    print(f"  (Groups: {', '.join(groups_pg_m)})")

    print("\nFinal Calculation:")
    print(f"Total number = {count_pg_2} (from point group '2') + {count_pg_m} (from point group 'm') = {total_count}")


if __name__ == "__main__":
    count_wallpaper_groups_by_point_group_order()