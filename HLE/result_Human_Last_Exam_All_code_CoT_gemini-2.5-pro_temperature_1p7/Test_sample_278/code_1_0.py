import math

def solve_area():
    """
    Calculates the area of the smallest square that outcircles the given parametric curve.
    """
    # Part 1: t in [0, pi], x = cos(t), y = sin(t)
    # x range is [-1, 1]
    # y range is [0, 1]
    x_min_1 = -1
    x_max_1 = 1
    y_min_1 = 0
    y_max_1 = 1

    # Part 2: t in (pi, 2pi), x = cos(t), y = 5*sin(t)
    # x range is [-1, 1]
    # y range is [-5, 0)
    x_min_2 = -1
    x_max_2 = 1
    y_min_2 = -5
    y_max_2 = 0

    # Overall range for the figure
    overall_x_min = min(x_min_1, x_min_2)
    overall_x_max = max(x_max_1, x_max_2)
    overall_y_min = min(y_min_1, y_min_2)
    overall_y_max = max(y_max_1, y_max_2)

    # Calculate the width and height of the bounding box
    width = overall_x_max - overall_x_min
    height = overall_y_max - overall_y_min

    # The side of the smallest square that outcircles the figure is the max of width and height
    side_length = max(width, height)

    # Calculate the area of the square
    area = side_length ** 2

    print(f"The total width of the figure is: {overall_x_max} - ({overall_x_min}) = {width}")
    print(f"The total height of the figure is: {overall_y_max} - ({overall_y_min}) = {height}")
    print(f"The side length of the smallest outcircling square is the maximum of the width and height: {side_length}")
    print(f"The area of the square is {int(side_length)} * {int(side_length)} = {int(area)}")

solve_area()