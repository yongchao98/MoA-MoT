import math

def solve_area():
    """
    Calculates the area of the smallest square that outcircles the given parametric curve.
    """
    # For t in [0, pi], x = cos(t), y = sin(t)
    # The range of cos(t) for t in [0, pi] is [-1, 1].
    # The range of sin(t) for t in [0, pi] is [0, 1].
    min_x_part1 = -1.0
    max_x_part1 = 1.0
    min_y_part1 = 0.0
    max_y_part1 = 1.0

    # For t in (pi, 2pi), x = cos(t), y = 5*sin(t)
    # The range of cos(t) for t in (pi, 2pi) is [-1, 1].
    # The range of sin(t) for t in (pi, 2pi) is [-1, 0), so y is in [-5, 0).
    min_x_part2 = -1.0
    max_x_part2 = 1.0
    min_y_part2 = -5.0
    max_y_part2 = 0.0

    # Overall range for the entire figure
    min_x = min(min_x_part1, min_x_part2)
    max_x = max(max_x_part1, max_x_part2)
    min_y = min(min_y_part1, min_y_part2)
    max_y = max(max_y_part1, max_y_part2)

    # Calculate the width and height of the figure's bounding box
    width = max_x - min_x
    height = max_y - min_y

    # The side of the smallest enclosing square is the maximum of the width and height
    side_length = max(width, height)

    # The area of the square
    area = side_length ** 2

    print(f"The range for x is [{min_x}, {max_x}].")
    print(f"The range for y is [{min_y}, {max_y}].")
    print(f"The width of the figure is {max_x} - ({min_x}) = {width}.")
    print(f"The height of the figure is {max_y} - ({min_y}) = {height}.")
    print(f"The side length of the smallest enclosing square is max({width}, {height}) = {side_length}.")
    print("The area of the square is side * side.")
    print(f"Final calculation: {int(side_length)} * {int(side_length)} = {int(area)}")

solve_area()