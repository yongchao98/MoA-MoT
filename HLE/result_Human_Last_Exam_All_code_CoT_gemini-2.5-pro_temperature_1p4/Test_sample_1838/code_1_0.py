def find_central_pixel(x1, y1, x2, y2):
    """
    Finds the central pixel of a line from (x1, y1) to (x2, y2)
    using Bresenham's line algorithm.
    This implementation is simplified for lines with a slope between 0 and 1.
    """
    
    # Calculate the change in x and y
    dx = x2 - x1
    dy = y2 - y1

    # Initialize the current point and the list of pixels
    x, y = x1, y1
    pixels = []

    # Initialize the decision parameter for Bresenham's algorithm
    # P = 2*dy - dx
    p = 2 * dy - dx
    
    print("Generating pixels using Bresenham's algorithm:")
    # Loop from the start x-coordinate to the end x-coordinate
    while x <= x2:
        # Add the current pixel to our list
        pixels.append((x, y))
        print(f"Pixel: ({x}, {y})")

        # Update x to move to the next column
        x += 1
        
        # Update the decision parameter and potentially y
        if p < 0:
            # P = P + 2*dy
            p = p + 2 * dy
        else:
            # P = P + 2*dy - 2*dx
            p = p + 2 * dy - 2 * dx
            y += 1
            
    # Find the middle index of the list of pixels
    center_index = (len(pixels) - 1) // 2
    
    # Get the central pixel using the index
    central_pixel = pixels[center_index]
    
    print("\n---")
    print(f"A total of {len(pixels)} pixels were generated.")
    print(f"The central pixel is the one at index {center_index}.")
    print(f"The coordinates of the central pixel are: {central_pixel[0]}, {central_pixel[1]}")

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

find_central_pixel(start_point[0], start_point[1], end_point[0], end_point[1])