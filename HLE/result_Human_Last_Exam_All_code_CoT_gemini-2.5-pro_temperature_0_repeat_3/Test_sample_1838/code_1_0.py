def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the pixels for a line using Bresenham's algorithm
    and returns the central pixel.
    This implementation is for slopes between 0 and 1.
    """
    pixels = []
    dx = x2 - x1
    dy = y2 - y1
    
    # Decision parameter
    D = 2 * dy - dx
    y = y1

    for x in range(x1, x2 + 1):
        pixels.append((x, y))
        # Check if the midpoint is above or below the line
        if D > 0:
            y += 1
            D += 2 * dy - 2 * dx
        else:
            D += 2 * dy
            
    # Find the central pixel
    if not pixels:
        print("No pixels were generated.")
        return

    middle_index = len(pixels) // 2
    central_pixel = pixels[middle_index]
    
    print(f"The line from ({x1},{y1}) to ({x2},{y2}) consists of {len(pixels)} pixels:")
    print(pixels)
    print("\nThe central pixel is the one at index " + str(middle_index) + " in the list.")
    print(f"Coordinates of the central pixel: {central_pixel}")

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

find_central_pixel_bresenham(start_point[0], start_point[1], end_point[0], end_point[1])

# The final answer is the coordinate tuple of the central pixel.
# The list of pixels is: [(1, 1), (2, 1), (3, 2), (4, 2), (5, 3), (6, 3), (7, 3), (8, 4), (9, 4), (10, 5), (11, 5)]
# There are 11 pixels, so the middle index is 11 // 2 = 5.
# The pixel at index 5 is (6, 3).