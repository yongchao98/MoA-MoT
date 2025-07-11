def find_central_pixel():
    """
    Calculates and prints the central pixel of a line using Bresenham's algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # The change in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Store the calculated pixels of the line
    pixels = []

    # Initial decision parameter
    # This version is for lines with slope between 0 and 1
    D = 2 * dy - dx
    y = y1

    # Iterate from start to end x-coordinate
    for x in range(x1, x2 + 1):
        pixels.append((x, y))
        
        # Check the decision parameter to decide the next y
        if D > 0:
            y += 1
            D += (2 * dy - 2 * dx)
        else:
            D += (2 * dy)

    # Calculate the index of the central pixel
    # For a list of 11 items, the middle index (0-based) is 5
    central_index = (len(pixels) - 1) // 2
    central_pixel = pixels[central_index]

    print(f"The line starts at ({x1}, {y1}) and ends at ({x2}, {y2}).")
    print(f"The pixels on the line are: {pixels}")
    print(f"The central pixel is the {central_index + 1}th pixel in the list.")
    print(f"The coordinates of the central pixel are: {central_pixel}")

find_central_pixel()
<<<(6, 3)>>>