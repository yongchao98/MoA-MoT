import math

def solve_area():
    """
    Calculates the area of the smallest square that outcircles the given parametric curve.
    """
    # For t in [0, pi]: x = cos(t), y = sin(t)
    # Range of x is [cos(pi), cos(0)] = [-1, 1]
    # Range of y is [sin(0 or pi), sin(pi/2)] = [0, 1]
    x_min_1 = -1
    x_max_1 = 1
    y_min_1 = 0
    y_max_1 = 1

    # For t in (pi, 2*pi): x = cos(t), y = 5*sin(t)
    # Range of x is [cos(pi), cos(2*pi)] = [-1, 1]
    # Range of y is [5*sin(3*pi/2), 5*sin(pi or 2*pi)] = [-5, 0]
    x_min_2 = -1
    x_max_2 = 1
    y_min_2 = -5
    y_max_2 = 0

    # Overall min and max values
    x_min = min(x_min_1, x_min_2)
    x_max = max(x_max_1, x_max_2)
    y_min = min(y_min_1, y_min_2)
    y_max = max(y_max_1, y_max_2)
    
    print(f"The figure spans from x={x_min} to x={x_max}.")
    print(f"The figure spans from y={y_min} to y={y_max}.")

    # Calculate the width and height of the bounding box
    width = x_max - x_min
    height = y_max - y_min

    print(f"\nThe width of the figure is {x_max} - ({x_min}) = {width}.")
    print(f"The height of the figure is {y_max} - ({y_min}) = {height}.")

    # The side of the smallest enclosing square is the max of width and height
    side_length = max(width, height)
    print(f"\nThe side length of the smallest enclosing square is max({width}, {height}) = {side_length}.")

    # Calculate the area of the square
    area = side_length * side_length
    print(f"\nThe area of the square is {side_length} * {side_length} = {area}.")

solve_area()