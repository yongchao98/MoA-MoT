import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the maximum and minimum value nodes
    in the described equilateral triangle problem.

    The problem can be solved by analyzing the properties of the values on the nodes.
    The rhombus condition implies that the value at any node is a linear function
    of its position. A linear function on a triangle attains its extrema at the
    vertices. Therefore, the max/min values are simply max/min of (a, b, c).
    The distance then depends on whether the vertices corresponding to the max
    and min values are distinct points, a point and a line, or the same set of points.
    """
    try:
        # Note: The problem statement gives 'n', but as shown in the derivation,
        # it is not needed for the final answer. We ask for a, b, c.
        print("Please enter the values a, b, and c for the vertices, separated by spaces.")
        a, b, c = map(float, input("Enter a b c: ").split())
    except ValueError:
        print("Invalid input. Please enter three numbers separated by spaces.")
        return

    # Use a set to find the number of unique values among a, b, c
    unique_values = set([a, b, c])
    num_unique = len(unique_values)

    if num_unique == 1:
        # Case: a = b = c
        # All nodes have the same value. The set of max points and min points is the same.
        print("\nAll three vertex values are equal.")
        print("The maximum and minimum values are the same, and all nodes have this value.")
        print("The shortest distance between the location of the max and min is 0.")
        print("Final Answer Equation: Distance = 0")

    elif num_unique == 2:
        # Case: Two values are equal, one is different. e.g., a = b != c
        # One extremum is at a vertex, the other is along the opposite side.
        # The distance is the altitude of the equilateral triangle.
        distance = math.sqrt(3) / 2
        print("\nExactly two of the vertex values are equal.")
        print("One extremum value is at a single vertex, and the other is along the opposite edge.")
        print("The shortest distance is the altitude of the triangle with side length 1.")
        print(f"Final Answer Equation: Distance = sqrt(3) / 2 = {distance}")

    else: # num_unique == 3
        # Case: a, b, c are all distinct.
        # The max is at one vertex and the min is at another.
        # The distance is the side length of the triangle.
        print("\nThe three vertex values are all distinct.")
        print("The maximum value is at one vertex and the minimum value is at another vertex.")
        print("The shortest distance is the side length between two vertices.")
        print("Final Answer Equation: Distance = 1")

if __name__ == '__main__':
    solve_triangle_distance()