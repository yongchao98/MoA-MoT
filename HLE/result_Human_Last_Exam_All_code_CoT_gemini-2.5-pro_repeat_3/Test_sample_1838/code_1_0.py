def find_central_pixel_bresenham(x0, y0, x1, y1):
    """
    Calculates the central pixel of a line using Bresenham's algorithm.
    This implementation is simplified for the specific case where dx > dy > 0.
    """
    # Calculate deltas
    dx = x1 - x0
    dy = y1 - y0

    # Initialize variables
    x, y = x0, y0
    p = 2 * dy - dx  # Initial decision parameter
    pixels = []

    # Generate all pixels on the line
    # The loop will run dx + 1 times
    while x <= x1:
        pixels.append((x, y))
        if p >= 0:
            y += 1
            p += 2 * dy - 2 * dx
        else:
            p += 2 * dy
        x += 1

    # Find the central pixel from the list
    middle_index = len(pixels) // 2
    central_pixel = pixels[middle_index]

    # Output the result
    # The final 'equation' is the coordinate pair (x, y)
    print(f"The coordinates of the central pixel are:")
    print(f"x = {central_pixel[0]}")
    print(f"y = {central_pixel[1]}")


# Define the start and end points of the line
start_x, start_y = 1, 1
end_x, end_y = 11, 5

# Run the function to find and print the central pixel
find_central_pixel_bresenham(start_x, start_y, end_x, end_y)