def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Finds the central pixel of a line using Bresenham's algorithm.
    This implementation is for lines with a positive slope less than 1.
    """
    # Calculate deltas
    dx = x2 - x1
    dy = y2 - y1

    # Initialize decision parameter
    # This is for a shallow slope (dx > dy)
    p = 2 * dy - dx

    # List to store the coordinates of all pixels on the line
    line_pixels = []

    # Initialize current coordinates
    x, y = x1, y1

    # Add the starting point to the list
    line_pixels.append((x, y))

    # Iterate from start x to end x
    while x < x2:
        x += 1
        if p < 0:
            # Move to the next pixel without changing y
            p = p + 2 * dy
        else:
            # Move to the next pixel and increment y
            y += 1
            p = p + 2 * dy - 2 * dx
        
        line_pixels.append((x, y))

    # The list of pixels is now generated. Find the central one.
    num_pixels = len(line_pixels)
    central_index = num_pixels // 2
    central_pixel = line_pixels[central_index]

    # Print the result
    print(f"All pixels on the line: {line_pixels}")
    print(f"Total number of pixels: {num_pixels}")
    print(f"The central pixel is at index {central_index}, which is the {central_index + 1}th pixel.")
    # The final print is formatted as per instructions.
    # The prompt asks to "output each number in the final equation", which means printing the coordinates.
    print(f"The central pixel's coordinates are: ({central_pixel[0]}, {central_pixel[1]})")

# Define the start and end points of the line
x1, y1 = 1, 1
x2, y2 = 11, 5

# Run the function to find and print the central pixel
find_central_pixel_bresenham(x1, y1, x2, y2)