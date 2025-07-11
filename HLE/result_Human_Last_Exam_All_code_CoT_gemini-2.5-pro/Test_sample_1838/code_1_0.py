def find_central_pixel():
    """
    Calculates the coordinates of the central pixel on a line from (1,1) to (11,5)
    using Bresenham's line algorithm.
    """
    # Define the start and end points of the line
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # Calculate the change in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Store all pixels of the line in a list
    pixels = []
    
    # Initialize current coordinates and the decision parameter
    x, y = x1, y1
    p = 2 * dy - dx
    
    # Pre-calculate constants for the loop
    const1 = 2 * dy
    const2 = 2 * (dy - dx)

    # Add the starting pixel
    pixels.append((x, y))

    # Iterate from the start to the end point
    while x < x2:
        x += 1
        if p < 0:
            p += const1
            # y does not change
        else:
            p += const2
            y += 1
        pixels.append((x, y))

    # Find the central pixel in the list of pixels
    # We use integer division to get the middle index
    mid_index = len(pixels) // 2
    central_pixel = pixels[mid_index]

    # Print the coordinates of the central pixel
    # This fulfills the requirement to "output each number in the final equation"
    print(f"The coordinates of the central pixel are: ({central_pixel[0]}, {central_pixel[1]})")

find_central_pixel()