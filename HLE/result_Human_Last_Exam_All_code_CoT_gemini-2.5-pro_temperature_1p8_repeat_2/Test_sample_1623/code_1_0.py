def solve_thurston_bennequin():
    """
    Calculates the Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # 1-based coordinates from the problem description
    o_coords_1based = [(1,1), (2,2), (3,3), (4,4), (5,5)]
    x_coords_1based = [(1,4), (2,5), (3,1), (4,2), (5,3)]

    # Convert to 0-based coordinates for easier processing
    o_coords = [(x - 1, y - 1) for x, y in o_coords_1based]
    x_coords = [(x - 1, y - 1) for x, y in x_coords_1based]
    
    all_coords = o_coords + x_coords

    nw_count = 0
    se_count = 0

    # Iterate through all pairs of (o, x)
    for i_o, j_o in o_coords:
        for i_x, j_x in x_coords:
            
            # According to the formula, we only need to consider pairs where the o is to the left of x
            if not (i_o < i_x):
                continue
            
            # Determine rectangle boundaries
            x_min, x_max = i_o, i_x
            y_min, y_max = min(j_o, j_x), max(j_o, j_x)

            # Check if the rectangle is empty
            is_empty = True
            for i_p, j_p in all_coords:
                # Skip the points defining the rectangle itself
                if (i_p, j_p) == (i_o, j_o) or (i_p, j_p) == (i_x, j_x):
                    continue
                
                # Check if the point is strictly inside the rectangle
                if x_min < i_p < x_max and y_min < j_p < y_max:
                    is_empty = False
                    break
            
            if is_empty:
                # Check for North-West configuration
                if j_o > j_x:
                    nw_count += 1
                # Check for South-East configuration
                elif j_o < j_x:
                    se_count += 1

    # Calculate the Thurston-Bennequin number
    tb_number = nw_count - se_count

    # Print the final equation
    print(f"Number of North-West corners (NW) = {nw_count}")
    print(f"Number of South-East corners (SE) = {se_count}")
    print(f"The maximal Thurston-Bennequin number = NW - SE = {nw_count} - {se_count} = {tb_number}")

solve_thurston_bennequin()