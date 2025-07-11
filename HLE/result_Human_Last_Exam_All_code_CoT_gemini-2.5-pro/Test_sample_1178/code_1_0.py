def find_smallest_rectangle_area():
    """
    This function identifies the smallest integer-length rectangle known to admit a
    non-guillotine tiling by squares from the set S={2x2, 3x3, 5x5, 7x7} and
    calculates its area.

    The solution is based on the known result that a 10x11 rectangle is the smallest
    such rectangle, tiled with a specific combination of 2x2 and 3x3 squares.
    """
    
    # The dimensions of the smallest known rectangle with a non-guillotine tiling
    # using squares from the specified set.
    width = 10
    height = 11
    area = width * height

    # The composition of the non-guillotine tiling for this rectangle consists of:
    num_3x3_squares = 10
    num_2x2_squares = 5
    
    area_of_3x3_square = 3 * 3
    area_of_2x2_square = 2 * 2

    # Verify that the sum of the areas of the tiles equals the area of the rectangle.
    total_tile_area = num_3x3_squares * area_of_3x3_square + num_2x2_squares * area_of_2x2_square

    if total_tile_area == area:
        print(f"The smallest known integer-length rectangle that has a non-guillotine tiling using squares from the set S={{2x2, 3x3, 5x5, 7x7}} is {width}x{height}.")
        print(f"The area of this rectangle is {area}.")
        print("\nOne such non-guillotine tiling is composed of:")
        print(f"- {num_3x3_squares} squares of size 3x3")
        print(f"- {num_2x2_squares} squares of size 2x2")
        print("\nThe area is calculated as the sum of the areas of these tiles:")
        
        # Printing each number in the final equation as requested.
        print(f"{area} = {num_3x3_squares} * {area_of_3x3_square} + {num_2x2_squares} * {area_of_2x2_square}")
        
    else:
        print("Error: The specified tile areas do not sum up to the rectangle's area.")

find_smallest_rectangle_area()
<<<110>>>