import math

def find_central_pixel_bresenham():
    """
    Finds the central pixel of a line using Bresenham's algorithm.
    """
    # Define the start and end points of the line
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # Lists to store the coordinates
    points = []

    # Calculate differences and initialize variables
    x, y = x1, y1
    dx = x2 - x1
    dy = y2 - y1

    # This implementation is for the gentle slope case (0 <= slope <= 1)
    # which applies here since dx > dy.
    
    # Initialize the decision parameter
    p = 2 * dy - dx

    # Loop through each x position to find the corresponding y
    while x <= x2:
        # Add the current pixel to our list
        points.append((x, y))

        # Update x for the next iteration
        x += 1
        
        # Update the decision parameter and y based on its value
        if p < 0:
            p = p + (2 * dy)
        else:
            p = p + (2 * dy) - (2 * dx)
            y += 1

    # Find the central pixel from the list of points
    # Using integer division to find the middle index
    mid_index = len(points) // 2
    central_pixel = points[mid_index]
    
    cx, cy = central_pixel
    
    # Print the final result, showing each number in the coordinate
    print(f"The coordinates of the central pixel are:")
    print(f"x = {cx}")
    print(f"y = {cy}")

# Run the function to find and print the central pixel
find_central_pixel_bresenham()