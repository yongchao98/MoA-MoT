import math

def find_shortest_distance():
    """
    Calculates the shortest distance between the max and min points based on the problem's logic.

    The problem reduces to finding the distance between vertices, or between a vertex and an edge
    of the equilateral triangle. The most non-trivial case, which gives a single numerical answer,
    is the distance from a vertex to the opposite side (the altitude).
    """
    side_length = 1.0
    half_side = side_length / 2.0

    # The shortest distance corresponds to the altitude of the equilateral triangle.
    # We calculate it using the Pythagorean theorem: altitude^2 + (side/2)^2 = side^2
    # The numbers in the final equation are the side length and half the side length.
    altitude = math.sqrt(side_length**2 - half_side**2)

    print(f"The final equation for the distance is: distance = sqrt({side_length}^2 - {half_side}^2)")
    print(f"The shortest distance is: {altitude}")

find_shortest_distance()