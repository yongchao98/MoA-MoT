def find_central_pixel():
    """
    Calculates the central pixel of a line using Bresenham's algorithm
    and prints the result.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # Calculate differences and initialize variables for Bresenham's algorithm.
    # This implementation is simplified for the case where the slope is between 0 and 1.
    dx = x2 - x1
    dy = y2 - y1
    p = 2 * dy - dx
    x, y = x1, y1

    pixels = []

    # The 'equation' for the line is the set of pixels it's composed of.
    # We will generate and store each pixel coordinate.
    while x <= x2:
        pixels.append((x, y))
        if p < 0:
            # Move to the next x, keep y the same. Update decision parameter.
            p = p + 2 * dy
        else:
            # Move to the next x and next y. Update decision parameter.
            p = p + 2 * dy - 2 * dx
            y = y + 1
        x = x + 1

    # Find the central pixel from the list of generated pixels.
    # Using integer division `//` finds the middle index.
    middle_index = len(pixels) // 2
    central_pixel = pixels[middle_index]

    print("Applying Bresenham's algorithm for the line from (1,1) to (11,5).")
    print("The complete set of pixels (the 'equation' of the line) is:")
    
    # We output each point (pair of numbers) that forms the line.
    points_str = ", ".join(map(str, pixels))
    print(points_str)
    
    print(f"\nThere are a total of {len(pixels)} pixels on the line.")
    print(f"The central pixel is the one at index {middle_index}, which is the {middle_index + 1}th pixel in the sequence.")
    print(f"The coordinates of the central pixel are: {central_pixel[0]}, {central_pixel[1]}")


find_central_pixel()