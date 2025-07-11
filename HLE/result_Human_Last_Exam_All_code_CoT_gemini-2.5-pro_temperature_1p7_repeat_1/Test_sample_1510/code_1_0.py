import collections

def solve_wallpaper_groups():
    """
    Calculates the number of wallpaper groups with a point group of order 2.
    """
    # In crystallography, there are 17 wallpaper groups (2D plane symmetry groups).
    # Each wallpaper group is associated with one of 10 possible 2D crystallographic point groups.
    # The 'order' of a point group is the number of symmetry elements it contains.

    # Data mapping the 17 wallpaper groups (IUC notation) to their point groups (Schönflies notation).
    wallpaper_to_point_group = {
        'p1': 'C1', 'p2': 'C2', 'pm': 'D1', 'pg': 'D1', 'cm': 'D1',
        'pmm': 'D2', 'pmg': 'D2', 'pgg': 'D2', 'cmm': 'D2', 'p3': 'C3',
        'p3m1': 'D3', 'p31m': 'D3', 'p4': 'C4', 'p4m': 'D4', 'p4g': 'D4',
        'p6': 'C6', 'p6m': 'D6'
    }

    # The order of each of the 10 2D point groups.
    # We are interested in point groups of order 2. These are C2 and D1.
    point_group_orders = {
        'C1': 1, 'C2': 2, 'D1': 2, 'C3': 3, 'C4': 4,
        'D2': 4, 'C6': 6, 'D3': 6, 'D4': 8, 'D6': 12
    }
    
    # We will count how many wallpaper groups belong to each point group of order 2.
    groups_by_point_group_order_2 = collections.defaultdict(list)
    for wg, pg in wallpaper_to_point_group.items():
        if point_group_orders.get(pg) == 2:
            groups_by_point_group_order_2[pg].append(wg)

    print("The point groups of order 2 are C2 and D1.")
    print("The wallpaper groups associated with these point groups are:")
    
    counts = []
    # Print the counts for each relevant point group.
    for pg in sorted(groups_by_point_group_order_2.keys()):
        count = len(groups_by_point_group_order_2[pg])
        counts.append(count)
        group_list_str = ", ".join(sorted(groups_by_point_group_order_2[pg]))
        print(f"- There are {count} groups with point group {pg}: {group_list_str}")

    # Display the final calculation as an equation.
    total = sum(counts)
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    
    print(f"\nThe total number of wallpaper groups with a point group of order 2 is:")
    print(f"{equation_str} = {total}")

solve_wallpaper_groups()