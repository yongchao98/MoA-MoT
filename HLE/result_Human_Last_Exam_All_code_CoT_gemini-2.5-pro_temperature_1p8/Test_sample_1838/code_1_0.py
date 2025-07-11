def find_central_pixel_bresenham():
    """
    Finds the central pixel of a line from (1, 1) to (11, 5)
    using Bresenham's line algorithm.
    """
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # List to store all the points on the line
    points = []

    # Calculate differences
    dx = x2 - x1
    dy = y2 - y1

    # Initialize current coordinates and decision parameter
    x, y = x1, y1
    p = 2 * dy - dx

    # Add the starting point
    points.append((x, y))

    # Loop to generate all points along the line
    # This loop will run dx times (10 times)
    while x < x2:
        x += 1
        if p < 0:
            # y remains the same, update decision parameter
            p = p + 2 * dy
        else:
            # y is incremented, update decision parameter
            y += 1
            p = p + 2 * dy - 2 * dx
        points.append((x, y))
    
    # Calculate the middle index
    mid_index = len(points) // 2
    
    # Get the central pixel coordinates
    central_pixel = points[mid_index]
    
    # As requested, output the numbers in the final answer
    print(f"The coordinates of the central pixel are: ({central_pixel[0]}, {central_pixel[1]})")

# Execute the function to find and print the result
find_central_pixel_bresenham()