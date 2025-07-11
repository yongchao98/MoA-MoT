def find_central_pixel():
    """
    Calculates the central pixel of a line using Bresenham's algorithm.
    """
    # Define the start and end points of the line
    x1, y1 = 1, 1
    x2, y2 = 11, 5

    # This list will store all the pixel coordinates on the line
    points = []

    # Calculate the change in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Initialize the current coordinates and the decision parameter 'p'
    # This specific formula for 'p' is for lines with a slope between 0 and 1
    x, y = x1, y1
    p = 2 * dy - dx

    # Add the starting point to our list
    points.append((x, y))

    # Loop for dx steps to find all the pixels along the line
    while x < x2:
        x += 1
        if p < 0:
            # Keep y the same
            p = p + (2 * dy)
        else:
            # Increment y
            y += 1
            p = p + (2 * dy) - (2 * dx)
        
        points.append((x, y))

    # Determine the index of the central pixel
    # For a list of 11 items (0-10), the middle index is (11-1)//2 = 5
    central_index = (len(points) - 1) // 2
    central_pixel = points[central_index]

    # The prompt asks to "output each number in the final equation".
    # We will print the x and y coordinates of the central pixel to satisfy this.
    print(f"The coordinates of the central pixel are (x, y).")
    print(f"x = {central_pixel[0]}")
    print(f"y = {central_pixel[1]}")

# Execute the function to find and print the result
find_central_pixel()