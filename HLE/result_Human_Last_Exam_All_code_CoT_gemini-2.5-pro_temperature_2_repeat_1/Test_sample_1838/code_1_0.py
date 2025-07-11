def find_central_pixel():
    """
    Calculates and prints the central pixel of a line using Bresenham's algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    pixels = []
    
    # Set up initial variables for Bresenham's algorithm
    x, y = x1, y1
    dx = x2 - x1
    dy = y2 - y1
    
    # This specific case has a slope between 0 and 1 (dx > dy > 0)
    # Initialize the decision parameter
    p = 2 * dy - dx

    # Iterate from start to end point to generate all pixels
    # We loop dx + 1 times to include the start and end points
    for _ in range(dx + 1):
        pixels.append((x, y))
        
        # Decide whether to increment y
        if p >= 0:
            y += 1
            p += 2 * dy - 2 * dx
        else:
            p += 2 * dy
        
        # Always increment x in this case
        x += 1

    # Find the central pixel from the list
    # The middle index is the total number of items // 2
    central_index = len(pixels) // 2
    central_pixel = pixels[central_index]

    print(f"The line from ({x1},{y1}) to ({x2},{y2}) contains {len(pixels)} pixels.")
    #Uncomment the line below to see all the pixels on the line
    #print(f"Pixels on the line: {pixels}")
    print(f"The central pixel is at coordinate: ({central_pixel[0]}, {central_pixel[1]})")

find_central_pixel()