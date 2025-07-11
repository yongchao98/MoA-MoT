def find_central_pixel_bresenham():
    """
    Finds the central pixel of a line using Bresenham's algorithm.
    The line is hardcoded from (1, 1) to (11, 5).
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    pixels = []
    x, y = x1, y1
    dx = x2 - x1
    dy = y2 - y1

    # This is a simplified version for slopes between 0 and 1.
    # Decision parameter
    p = 2 * dy - dx

    # Add the starting point
    pixels.append((x, y))

    # Iterate from x1 to x2-1
    while x < x2:
        x += 1
        if p < 0:
            p = p + 2 * dy
            # y does not change
        else:
            p = p + 2 * dy - 2 * dx
            y += 1
        pixels.append((x, y))
    
    # Find the middle index of the list of pixels.
    # For a list with an odd number of elements, this gives the exact center.
    middle_index = len(pixels) // 2
    central_pixel = pixels[middle_index]

    # Output the coordinates of the central pixel
    print(f"The coordinates of the central pixel are ({central_pixel[0]}, {central_pixel[1]})")

find_central_pixel_bresenham()
<<< (6, 3) >>>