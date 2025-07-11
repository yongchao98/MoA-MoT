def calculate_intersections(ux, uy):
    """
    Calculates the intersection of a unit square U(ux, uy) with
    two polygons: a "black" and a "white" set of squares forming a
    checkerboard on a 4x4 grid.

    Args:
        ux: The x-coordinate of the bottom-left corner of the unit square.
        uy: The y-coordinate of the bottom-left corner of the unit square.
    """
    print(f"Analyzing the unit square U with bottom-left corner at ({ux}, {uy}).\n")

    area_black = 0.0
    area_white = 0.0
    
    # The checkerboard polygons are unions of the 16 unit squares (cells).
    # We iterate through each cell, calculate the intersection with U,
    # and add it to the corresponding polygon's total area.
    print("Decomposition is a 2-polygon checkerboard.")
    print("Black polygon is the union of all cells (i,j) where i+j is even.")
    print("White polygon is the union of all cells (i,j) where i+j is odd.\n")

    print("Calculating intersection with the Black polygon:")
    # Iterate through the 16 cells of the 4x4 grid
    for i in range(4):
        for j in range(4):
            # Check color of the cell (i, j)
            is_black = (i + j) % 2 == 0
            
            # Calculate the intersection of the unit square U and the cell
            x_overlap = max(0, min(ux + 1, i + 1) - max(ux, i))
            y_overlap = max(0, min(uy + 1, j + 1) - max(uy, j))
            intersection_area = x_overlap * y_overlap
            
            if intersection_area > 0:
                if is_black:
                    print(f"  Intersection with cell ({i},{j}) (Black): {x_overlap:.2f} * {y_overlap:.2f} = {intersection_area:.4f}")
                    area_black += intersection_area
                else: # is_white
                    pass # Will calculate in the next block

    print(f"Total Area(U intersect Black Polygon) = {area_black:.4f}\n")

    print("Calculating intersection with the White polygon:")
    for i in range(4):
        for j in range(4):
            is_black = (i + j) % 2 == 0
            x_overlap = max(0, min(ux + 1, i + 1) - max(ux, i))
            y_overlap = max(0, min(uy + 1, j + 1) - max(uy, j))
            intersection_area = x_overlap * y_overlap
            
            if intersection_area > 0:
                if not is_black:
                    print(f"  Intersection with cell ({i},{j}) (White): {x_overlap:.2f} * {y_overlap:.2f} = {intersection_area:.4f}")
                    area_white += intersection_area
    
    print(f"Total Area(U intersect White Polygon) = {area_white:.4f}\n")

    total_area = area_black + area_white
    max_intersection = max(area_black, area_white)

    print("--- Summary ---")
    print(f"Total intersection area (sum) = {area_black:.4f} + {area_white:.4f} = {total_area:.4f} (should be 1.0)")
    print(f"Maximum intersection area = max({area_black:.4f}, {area_white:.4f}) = {max_intersection:.4f}")
    print("\nFor any unit square, this maximum is always >= 0.5.")
    print("This principle can be extended to a 16-polygon case, making the largest possible r = 0.5.")


# Example: A "worst-case" unit square for the checkerboard, centered on a grid line.
# Any U where a corner is at (i+0.5, j) or (i, j+0.5) gives a 0.5/0.5 split.
calculate_intersections(ux=0.5, uy=1.0)

<<<1/2>>>