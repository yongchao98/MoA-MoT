import math

def solve_area():
    """
    Calculates the area of the smallest square enclosing the parametric curve.
    """
    # Step 1: Determine the extents of the curve.
    # The x-coordinate is x = cos(t). For t in [0, 2pi), its range is [-1, 1].
    x_min = -1.0
    x_max = 1.0

    # The y-coordinate has two parts.
    # For t in [0, pi], y = 1 * sin(t). The range of sin(t) is [0, 1].
    y_range1_min = 0.0
    y_range1_max = 1.0

    # For t in (pi, 2pi), y = 5 * sin(t). The range of sin(t) is [-1, 0).
    y_range2_min = 5.0 * -1.0
    y_range2_max = 0.0

    # The total range for y is the union of the two ranges.
    y_min = min(y_range1_min, y_range2_min)
    y_max = max(y_range1_max, y_range2_max)

    print("Step 1: Determine the extents of the curve.")
    print(f"The range for x is [{int(x_min)}, {int(x_max)}].")
    print(f"The range for y is [{int(y_min)}, {int(y_max)}].")
    print("-" * 30)

    # Step 2: Calculate the dimensions of the bounding box.
    width = x_max - x_min
    height = y_max - y_min
    
    print("Step 2: Calculate the dimensions of the bounding box.")
    print(f"Width = {int(x_max)} - ({int(x_min)}) = {int(width)}")
    print(f"Height = {int(y_max)} - ({int(y_min)}) = {int(height)}")
    print("-" * 30)

    # Step 3: Find the side length of the smallest enclosing square.
    side_length = max(width, height)

    print("Step 3: Find the side length of the smallest enclosing square.")
    print(f"Side length = max({int(width)}, {int(height)}) = {int(side_length)}")
    print("-" * 30)

    # Step 4: Calculate the area of the square.
    area = side_length ** 2

    print("Step 4: Calculate the final area.")
    print(f"Area = {int(side_length)}^2 = {int(area)}")

solve_area()