def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the coordinates of the central pixel on a line using Bresenham's algorithm.
    """
    # Calculate deltas
    dx = x2 - x1
    dy = y2 - y1

    # Check for steep lines (slope > 1) or shallow lines (slope <= 1)
    # This specific problem has a shallow slope (dy/dx = 4/10 < 1)
    # The code below is simplified for this case where 0 < slope <= 1
    
    # List to store all pixels on the line
    pixels = []
    
    # Initialize current coordinates and decision parameter
    x, y = x1, y1
    p = 2 * dy - dx
    
    # Add the starting pixel
    pixels.append((x, y))
    
    # Iterate from start x to end x
    while x < x2:
        x += 1
        if p < 0:
            # y does not change
            p = p + 2 * dy
        else:
            # y increments
            y += 1
            p = p + 2 * dy - 2 * dx
        pixels.append((x, y))
        
    # Find the central pixel from the list
    mid_index = len(pixels) // 2
    central_pixel = pixels[mid_index]
    
    # Print the final result
    print(f"The line is between ({x1},{y1}) and ({x2},{y2}).")
    print(f"The pixels on the line are: {pixels}")
    print(f"The total number of pixels is: {len(pixels)}.")
    print(f"The central pixel is the one at index {mid_index}.")
    print(f"The coordinates of the central pixel are: {central_pixel}")

# Define the start and end points of the line
start_x, start_y = 1, 1
end_x, end_y = 11, 5

find_central_pixel_bresenham(start_x, start_y, end_x, end_y)
<<<(6, 3)>>>