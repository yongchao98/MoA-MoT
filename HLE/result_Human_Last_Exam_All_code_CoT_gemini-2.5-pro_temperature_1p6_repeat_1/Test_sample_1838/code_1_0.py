def find_central_pixel_bresenham(x1, y1, x2, y2):
    """
    Calculates the pixels for a line using Bresenham's algorithm
    and returns the central pixel. This implementation is for slopes < 1.
    """
    pixels = []
    dx = x2 - x1
    dy = y2 - y1

    # Decision parameter
    p = 2 * dy - dx
    x, y = x1, y1

    # Loop through each x coordinate
    while x <= x2:
        pixels.append((x, y))
        x += 1
        if p < 0:
            p = p + 2 * dy
        else:
            p = p + 2 * dy - 2 * dx
            y += 1
    
    print("All pixels on the line:", pixels)
    
    # Find the central pixel
    mid_index = len(pixels) // 2
    central_pixel = pixels[mid_index]
    
    print("\nTotal number of pixels:", len(pixels))
    print("The central pixel is at index:", mid_index)
    
    # Output each number in the final result
    print("\nThe coordinates of the central pixel are:")
    print(f"x = {central_pixel[0]}")
    print(f"y = {central_pixel[1]}")

    return central_pixel

# Define the start and end points of the line
start_point = (1, 1)
end_point = (11, 5)

# Find and print the central pixel
central_pixel_coords = find_central_pixel_bresenham(start_point[0], start_point[1], end_point[0], end_point[1])

# The final answer in the requested format
# The output format for the final answer should be a tuple (x, y)
# For example: (6, 3)
final_answer_string = f"({central_pixel_coords[0]}, {central_pixel_coords[1]})"
# print(f"\n<<<{final_answer_string}>>>")