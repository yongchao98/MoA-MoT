def find_central_pixel(x1, y1, x2, y2):
    """
    Calculates the central pixel of a line using Bresenham's algorithm.
    This implementation is for slopes between 0 and 1.
    """
    # Calculate the differences in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Initialize the decision parameter
    # p = 2*dy - dx
    p = 2 * dy - dx

    # Initialize starting coordinates
    x, y = x1, y1

    # Store the pixels of the line
    pixels = []

    # Iterate from x1 to x2
    while x <= x2:
        # Add the current pixel to the list
        pixels.append((x, y))

        # Update x for the next iteration
        x += 1

        # Check the decision parameter to decide the next y
        if p < 0:
            # If p is negative, the next y is the same. Update p.
            p = p + 2 * dy
        else:
            # If p is non-negative, increment y. Update p.
            y += 1
            p = p + 2 * dy - 2 * dx

    # Find the central pixel from the list
    # The middle index is found by integer division of the list length
    middle_index = len(pixels) // 2
    central_pixel = pixels[middle_index]

    print(f"The line from ({x1},{y1}) to ({x2},{y2}) has the following pixels:")
    print(pixels)
    print(f"\nThe total number of pixels is {len(pixels)}.")
    print(f"The central pixel is at index {middle_index}.")
    print(f"The coordinates of the central pixel are: {central_pixel}")

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

# Find and print the central pixel
find_central_pixel(start_point[0], start_point[1], end_point[0], end_point[1])

<<<(6, 3)>>>