def find_central_pixel(x1, y1, x2, y2):
    """
    Calculates all pixels on a line using Bresenham's algorithm and finds the central one.
    This implementation is for slopes between 0 and 1.
    """
    dx = x2 - x1
    dy = y2 - y1

    # Decision parameter
    p = 2 * dy - dx
    
    # Store points on the line
    points = []
    
    x, y = x1, y1

    # Loop for dx steps
    for _ in range(dx + 1):
        points.append((x, y))
        
        if x == x2:
            break
            
        x += 1
        if p < 0:
            # y doesn't change
            p = p + 2 * dy
        else:
            # y increments
            y += 1
            p = p + 2 * dy - 2 * dx

    print("All pixels on the line using Bresenham's algorithm:")
    print(points)
    print("-" * 30)

    # Find the central pixel
    mid_index = len(points) // 2
    central_pixel = points[mid_index]

    print(f"Total number of pixels: {len(points)}")
    print(f"The central pixel is the {mid_index + 1}th pixel in the list.")
    print(f"The coordinates of the central pixel are: {central_pixel}")

# Given line coordinates
start_point = (1, 1)
end_point = (11, 5)

find_central_pixel(start_point[0], start_point[1], end_point[0], end_point[1])
<<<(6, 3)>>>