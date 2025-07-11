def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the pixels on a line using Bresenham's algorithm
    and finds the central pixel.
    This implementation is for slopes between 0 and 1.
    """
    # List to store all the points on the line
    points = []

    # Calculate the differences in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Initialize the current point and the decision parameter
    x, y = x1, y1
    p = 2 * dy - dx

    # Add the starting point
    points.append((x, y))
    
    # Iterate from start x to end x
    while x < x2:
        x += 1
        if p < 0:
            # The next point is (x, y)
            p = p + 2 * dy
        else:
            # The next point is (x, y+1)
            y += 1
            p = p + 2 * dy - 2 * dx
        points.append((x, y))

    # Find the central pixel
    mid_index = (len(points) - 1) // 2
    central_pixel = points[mid_index]

    # Print all the pixels on the line
    print("The pixels on the line are:")
    for point in points:
        print(point)
    
    # Print the final result
    print("\n-----------------------------------------")
    print(f"The line has {len(points)} pixels in total.")
    print(f"The central pixel is the {mid_index + 1}th pixel.")
    print(f"The coordinates of the central pixel are: ({central_pixel[0]}, {central_pixel[1]})")

# Define the start and end points of the line
start_x, start_y = 1, 1
end_x, end_y = 11, 5

# Run the function
find_central_pixel_bresenham(start_x, start_y, end_x, end_y)
