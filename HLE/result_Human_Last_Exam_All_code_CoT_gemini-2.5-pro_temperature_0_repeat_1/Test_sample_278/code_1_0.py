import math

def solve_area():
    """
    Calculates the area of the smallest square that outcircles the given parametric curve.
    """
    # Step 1: Analyze the curve for t in [0, pi]
    # x = cos(t), y = sin(t)
    # For x = cos(t), max is at t=0 (cos(0)=1), min is at t=pi (cos(pi)=-1)
    # For y = sin(t), max is at t=pi/2 (sin(pi/2)=1), min is at t=0,pi (sin(0)=sin(pi)=0)
    x_max_1 = 1
    x_min_1 = -1
    y_max_1 = 1
    y_min_1 = 0
    
    print("Analyzing the first part of the curve (t in [0, pi]):")
    print(f"  x = cos(t) ranges from {x_min_1} to {x_max_1}")
    print(f"  y = sin(t) ranges from {y_min_1} to {y_max_1}")
    print("-" * 20)

    # Step 2: Analyze the curve for t in (pi, 2*pi)
    # x = cos(t), y = 5*sin(t)
    # For x = cos(t), max is at t=2pi (cos(2pi)=1), min is at t=pi (cos(pi)=-1)
    # For y = 5*sin(t), max is at t=pi,2pi (5*sin(pi)=0), min is at t=3pi/2 (5*sin(3pi/2)=-5)
    x_max_2 = 1
    x_min_2 = -1
    y_max_2 = 0
    y_min_2 = -5

    print("Analyzing the second part of the curve (t in (pi, 2*pi)):")
    print(f"  x = cos(t) ranges from {x_min_2} to {x_max_2}")
    print(f"  y = 5*sin(t) ranges from {y_min_2} to {y_max_2}")
    print("-" * 20)

    # Step 3: Find the overall bounding box for the entire figure
    x_min = min(x_min_1, x_min_2)
    x_max = max(x_max_1, x_max_2)
    y_min = min(y_min_1, y_min_2)
    y_max = max(y_max_1, y_max_2)

    print("Finding the overall bounding box:")
    print(f"  Overall x-range: [{x_min}, {x_max}]")
    print(f"  Overall y-range: [{y_min}, {y_max}]")
    print("-" * 20)

    # Step 4: Calculate the width and height of the bounding box
    width = x_max - x_min
    height = y_max - y_min
    
    print("Calculating the dimensions of the bounding rectangle:")
    print(f"  Width = {x_max} - ({x_min}) = {width}")
    print(f"  Height = {y_max} - ({y_min}) = {height}")
    print("-" * 20)

    # Step 5: Determine the side length of the smallest enclosing square
    side_length = max(width, height)
    
    print("Determining the side length of the smallest enclosing square:")
    print(f"  Side length = max(width, height) = max({width}, {height}) = {side_length}")
    print("-" * 20)

    # Step 6: Calculate the area of the square
    area = side_length ** 2
    
    print("Calculating the final area:")
    print(f"  Area = side_length * side_length")
    print(f"  Area = {side_length} * {side_length} = {area}")

solve_area()