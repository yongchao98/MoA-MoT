def find_central_pixel(x1, y1, x2, y2):
    """
    Calculates the pixels on a line using Bresenham's algorithm and finds the central one.
    """
    # Create a list to store all the pixels of the line
    pixels = []
    
    # Calculate the differences in x and y
    dx = abs(x2 - x1)
    dy = -abs(y2 - y1)
    
    # Determine the direction of the step for x and y
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    
    # Initialize the error term
    err = dx + dy
    
    # Set the starting point
    x, y = x1, y1
    
    # Loop to generate all pixels on the line
    while True:
        pixels.append((x, y))
        if x == x2 and y == y2:
            break
        
        # Calculate the error for the next step
        e2 = 2 * err
        
        # Update x coordinate and error term
        if e2 >= dy:
            err += dy
            x += sx
            
        # Update y coordinate and error term
        if e2 <= dx:
            err += dx
            y += sy
            
    # Find the middle index of the list of pixels
    mid_index = len(pixels) // 2
    
    # Get the central pixel using the middle index
    central_pixel = pixels[mid_index]
    
    print(f"The line from ({x1}, {y1}) to ({x2}, {y2}) has {len(pixels)} pixels.")
    print(f"The central pixel is at index {mid_index} (0-indexed).")
    print(f"The coordinates of the central pixel are: {central_pixel}")

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

# Run the function and print the result
find_central_pixel(start_point[0], start_point[1], end_point[0], end_point[1])