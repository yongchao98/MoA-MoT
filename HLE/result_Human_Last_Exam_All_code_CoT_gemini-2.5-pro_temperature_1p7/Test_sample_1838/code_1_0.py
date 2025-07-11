def find_central_pixel_bresenham():
    """
    Calculates and prints the central pixel of a line from (1, 1) to (11, 5)
    using Bresenham's line algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    pixels = []
    
    # These variables are for the algorithm
    x, y = x1, y1
    dx = x2 - x1
    dy = y2 - y1
    
    # Decision parameter
    p = 2 * dy - dx

    # Store the first pixel
    pixels.append((x, y))

    # Iterate from x1 to x2 to generate all pixels
    while x < x2:
        x += 1
        if p < 0:
            p = p + 2 * dy
        else:
            y += 1
            p = p + 2 * dy - 2 * dx
        pixels.append((x, y))
    
    # Find the middle pixel in the list of generated pixels
    mid_index = len(pixels) // 2
    central_pixel = pixels[mid_index]
    
    # Output the result as requested
    center_x = central_pixel[0]
    center_y = central_pixel[1]
    
    print(f"The central pixel coordinate is ({center_x}, {center_y})")

find_central_pixel_bresenham()