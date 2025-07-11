def find_central_pixel():
    """
    Calculates and prints the central pixel of a line using Bresenham's algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # Store all pixels of the line
    pixels = []

    # Calculate deltas
    dx = x2 - x1
    dy = y2 - y1

    # Decision parameter, initialized for the first step
    p = 2 * dy - dx

    x, y = x1, y1

    # The loop runs as long as x is less than or equal to the end x-coordinate
    while x <= x2:
        pixels.append((x, y))
        if p < 0:
            # Next point is (x+1, y)
            p = p + 2 * dy
        else:
            # Next point is (x+1, y+1)
            p = p + 2 * dy - 2 * dx
            y += 1
        x += 1
    
    # The loop generates one extra point, so we remove it.
    # A more common implementation would loop `dx` times.
    # Let's adjust the loop to be more standard.
    
    # --- Corrected Bresenham's Implementation ---
    pixels_corrected = []
    x, y = x1, y1
    dx = x2 - x1
    dy = y2 - y1
    p = 2 * dy - dx
    
    # Add the starting point
    pixels_corrected.append((x, y))
    
    # Loop dx times to generate the remaining points
    for _ in range(dx):
        x += 1
        if p < 0:
            p = p + 2 * dy
        else:
            y += 1
            p = p + 2 * dy - 2 * dx
        pixels_corrected.append((x, y))

    # Find the middle index
    # For a list of N items, the middle index is N // 2
    middle_index = len(pixels_corrected) // 2
    central_pixel = pixels_corrected[middle_index]

    # Output the final result
    print(f"The central pixel is at coordinates: ({central_pixel[0]}, {central_pixel[1]})")

find_central_pixel()
<<<(6, 3)>>>