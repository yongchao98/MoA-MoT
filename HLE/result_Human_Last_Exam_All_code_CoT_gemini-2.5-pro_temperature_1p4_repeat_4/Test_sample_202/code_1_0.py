import math

def find_shortest_distance():
    """
    This function calculates the shortest distance between the nodes with the maximum
    and minimum values on the described triangular grid.

    You can modify the values for a, b, and c below to test different inputs.
    """
    # ---- User-configurable values ----
    # These are the real numbers placed on the vertices A, B, and C.
    a = 10
    b = 20
    c = 10
    # ----------------------------------

    print(f"Calculating shortest distance for vertex values: a={a}, b={b}, c={c}\n")

    # The number of unique values among a, b, and c determines the solution.
    num_unique_values = len(set([a, b, c]))

    if num_unique_values == 1:
        # Case: a = b = c. All points on the grid have the same value.
        # The shortest distance between a point with the max value and one with the min value is 0.
        print("Analysis: All three vertex values are equal.")
        print("The maximum and minimum values are the same and co-located.")
        distance = 0
        print(f"\nThe shortest distance is {distance}.")
        print("Final equation: distance = 0")
        print("The number in the equation is: 0")

    elif num_unique_values == 3:
        # Case: a, b, and c are all different.
        # The max is at one vertex, and the min is at another. The distance is the side length.
        print("Analysis: The three vertex values are distinct.")
        print("The maximum is at one vertex and the minimum is at another.")
        distance = 1
        print(f"\nThe shortest distance is the side length of the triangle, which is {distance}.")
        print("Final equation: distance = 1")
        print("The number in the equation is: 1")

    else:  # num_unique_values == 2
        # Case: Exactly two values are equal.
        # One extremum is at a vertex, the other is on the opposite edge. Distance is the altitude.
        print("Analysis: Exactly two of the vertex values are equal.")
        print("One extremum is at a single vertex, and the other is on the opposite edge.")
        distance = math.sqrt(3) / 2
        print(f"\nThe shortest distance is the triangle's altitude, sqrt(3)/2.")
        print(f"Final equation: distance = sqrt(3) / 2")
        print(f"The numbers in the equation's expression are 3 and 2.")
        print(f"Numerical value: {distance:.5f}")


if __name__ == '__main__':
    find_shortest_distance()