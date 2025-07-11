def find_central_pixel(x1, y1, x2, y2):
    """
    Calculates the coordinates of the central pixel on a line
    using Bresenham's algorithm.
    """
    # Calculate deltas and initial decision parameter
    dx = x2 - x1
    dy = y2 - y1
    
    # Check if the line is steep or gentle
    is_steep = abs(dy) > abs(dx)
    
    # Rotate the line if it's steep
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2
    
    # Swap endpoints if necessary
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True
        
    # Recalculate deltas after potential swapping
    dx = x2 - x1
    dy = y2 - y1
    
    # Initialize error term and y-step direction
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1
    
    # Start at the beginning of the line
    y = y1
    points = []
    
    # Iterate over x coordinates
    for x in range(x1, x2 + 1):
        # Determine the coordinate based on rotation
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        
        # Update error term
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx
            
    # Reverse the list if endpoints were swapped
    if swapped:
        points.reverse()
        
    # Find and print the central pixel
    if not points:
        print("No pixels generated.")
        return

    center_index = (len(points) - 1) // 2
    central_pixel = points[center_index]
    
    print(f"Line from ({points[0][0]},{points[0][1]}) to ({points[-1][0]},{points[-1][1]})")
    print(f"Total pixels on the line: {len(points)}")
    print(f"The central pixel is number {center_index + 1} in the sequence.")
    print(f"The coordinates of the central pixel are: {central_pixel[0]}, {central_pixel[1]}")

# Define the start and end points of the line
start_x, start_y = 1, 1
end_x, end_y = 11, 5

find_central_pixel(start_x, start_y, end_x, end_y)

<<<(6, 3)>>>