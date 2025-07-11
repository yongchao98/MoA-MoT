def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the central pixel of a line using Bresenham's algorithm.
    This implementation is specific to the case where x1 < x2 and the slope is between 0 and 1.
    """
    # Calculate the change in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Initialize the decision parameter.
    # This is the initial error term.
    p = 2 * dy - dx

    # Initialize starting coordinates
    x, y = x1, y1

    # Store all pixels on the line in a list
    pixels = []

    # Iterate from the start to the end of the line
    while x <= x2:
        # Add the current pixel to our list
        pixels.append((x, y))

        # Bresenham's algorithm logic
        if p < 0:
            # The next point is (x+1, y)
            p = p + 2 * dy
        else:
            # The next point is (x+1, y+1)
            p = p + 2 * dy - 2 * dx
            y = y + 1
        
        # Always increment x for this slope
        x = x + 1

    # Find the central pixel
    # The number of pixels is len(pixels). The middle index is len(pixels) // 2.
    central_index = len(pixels) // 2
    central_pixel = pixels[central_index]

    # Print the coordinates of the central pixel
    print(f"The line from ({x1},{y1}) to ({x2},{y2}) contains {len(pixels)} pixels.")
    print("The list of all pixels is:", pixels)
    print("\nThe central pixel is the one at index", central_index)
    
    # The request is to output each number in the final equation.
    # The final "equation" is the coordinate pair of the central pixel.
    final_x = central_pixel[0]
    final_y = central_pixel[1]
    
    print(f"\nThe coordinates of the central pixel are ({final_x}, {final_y})")


# Given start and end points
start_point = (1, 1)
end_point = (11, 5)

find_central_pixel_bresenham(start_point[0], start_point[1], end_point[0], end_point[1])