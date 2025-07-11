import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the maximum and minimum value nodes
    in a divided equilateral triangle.

    The problem states an equilateral triangle ABC with side length 1 has values a, b, c
    at its vertices. A key property implies the value at any node is a linear
    interpolation of the vertex values. This means the maximum and minimum values
    for the entire triangle must lie at the vertices A, B, or C.
    """

    # --- User-configurable inputs ---
    # Set the values for the numbers placed on vertices A, B, and C.
    # Note: The side length is 1, and 'n' (the number of divisions) does not affect the final distance.
    a = 5.0
    b = 5.0
    c = 2.0
    # --- End of inputs ---

    print(f"Given values at vertices A, B, C are: a={a}, b={b}, c={c}")

    # The locations of the max/min values depend on how many of a,b,c are unique.
    unique_values = len(set([a, b, c]))

    if unique_values == 1:
        # Case 1: a = b = c.
        # The value is constant everywhere. The max and min points are the same.
        # The shortest distance between the set of max points and min points is 0.
        print("The three vertex values are identical.")
        final_distance = 0.0
        print(f"The shortest distance is: {final_distance}")

    elif unique_values == 2:
        # Case 2: Two values are equal, one is different (e.g., a=b > c).
        # One extreme value is on an edge, the other is at the opposite vertex.
        # The shortest distance is the altitude of the equilateral triangle.
        side_length = 1.0
        altitude = math.sqrt(side_length**2 - (side_length/2)**2)
        print("Two of the three vertex values are identical.")
        print("The shortest distance is the altitude of the triangle.")
        # The final equation for the distance:
        print(f"Distance = sqrt({side_length}^2 - ({side_length}/2)^2) = {altitude}")
        final_distance = altitude

    else:  # unique_values == 3
        # Case 3: a, b, c are all distinct.
        # The max is at one vertex, and the min is at another.
        # The shortest distance is the side length between these two vertices.
        side_length = 1.0
        print("The three vertex values are distinct.")
        print("The shortest distance is the side length of the triangle.")
        final_distance = side_length
        # The final equation for the distance:
        print(f"Distance = {final_distance}")

solve_triangle_distance()