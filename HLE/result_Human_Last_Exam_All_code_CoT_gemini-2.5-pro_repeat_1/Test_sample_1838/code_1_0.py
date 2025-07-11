def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the pixels for a line using Bresenham's algorithm
    and prints the coordinates of the central pixel.
    """
    # Calculate the differences in x and y
    dx = x2 - x1
    dy = y2 - y1

    # For this specific case, slope is 0 < m < 1, so we iterate along x
    
    # Initialize the decision parameter
    p = 2 * dy - dx
    
    # Initialize starting coordinates
    x, y = x1, y1
    
    # Store all pixels on the line
    pixels = []
    
    # Loop to generate all points until the end point is reached
    while x <= x2:
        pixels.append((x, y))
        
        # Increment x for the next step
        x += 1
        
        # Update the decision parameter and y coordinate
        if p < 0:
            p = p + 2 * dy
        else:
            p = p + 2 * dy - 2 * dx
            y += 1
            
    # Find the middle index
    mid_index = len(pixels) // 2
    
    # Get the central pixel from the list
    central_pixel = pixels[mid_index]
    
    # The final equation is finding the central pixel's coordinates.
    # The following line outputs each number (coordinate) of that final result.
    print(f"The coordinates of the central pixel are ({central_pixel[0]}, {central_pixel[1]})")

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

# Run the function
find_central_pixel_bresenham(start_point[0], start_point[1], end_point[0], end_point[1])