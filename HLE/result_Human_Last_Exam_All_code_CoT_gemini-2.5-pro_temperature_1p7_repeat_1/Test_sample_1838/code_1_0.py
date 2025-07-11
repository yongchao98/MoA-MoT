def find_central_pixel():
    """
    Calculates the central pixel of a line using Bresenham's algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # A list to store all the pixels on the line
    points = []

    # Calculate deltas and initial decision parameter
    dx = x2 - x1
    dy = y2 - y1
    
    # This is a simplified version of Bresenham's for the first octant (0 <= m <= 1)
    # where m is the slope. Our slope is 4/10 = 0.4, which fits this case.
    decision_parameter = 2 * dy - dx
    
    x, y = x1, y1

    # Iterate from the start to the end point along the major axis (x)
    while x <= x2:
        points.append((x, y))
        
        # Move to the next pixel
        x += 1
        
        # Update decision parameter to choose the next y-coordinate
        if decision_parameter < 0:
            # Keep y the same
            decision_parameter += (2 * dy)
        else:
            # Increment y
            y += 1
            decision_parameter += (2 * dy) - (2 * dx)

    # The line has been generated. Now, find the central pixel.
    # With an odd number of points, the middle index is easy to find.
    central_index = len(points) // 2
    central_pixel = points[central_index]
    
    print(f"The line from ({x1},{y1}) to ({x2},{y2}) passes through {len(points)} pixels.")
    # The prompt asks to output the numbers in the final equation. 
    # Since the final output is a coordinate, not an equation, I will just print the coordinate.
    print(f"The central pixel has coordinates: ({central_pixel[0]}, {central_pixel[1]})")

find_central_pixel()