import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the maximum and minimum value nodes
    on the described triangular grid.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    # The side length of the large equilateral triangle is 1.
    side_length = 1.0

    # The solution depends on the number of unique values among a, b, and c.
    # We use a set to find the number of unique values.
    unique_values_count = len(set([a, b, c]))

    print(f"Analyzing for vertex values a={a}, b={b}, c={c}:")

    if unique_values_count == 1:
        # Case 1: a = b = c.
        # All nodes have the same value. The set of maximum points is the same as
        # the set of minimum points. The shortest distance is 0.
        distance = 0.0
        print("All three vertex values are equal.")
        print(f"The shortest distance is: {distance}")

    elif unique_values_count == 2:
        # Case 2: Exactly two values are equal (e.g., a=b > c or a > b=c).
        # One extremum occurs at a single vertex, and the other occurs along the
        # entire opposite edge. The shortest distance is the altitude of the
        # equilateral triangle.
        # Altitude h = (sqrt(3)/2) * side_length
        distance = (math.sqrt(3) / 2) * side_length
        print("Exactly two of the vertex values are equal.")
        print(f"The shortest distance is the triangle's altitude.")
        print(f"Equation: sqrt(3)/2 * {side_length}")
        print(f"Result: {distance}")

    else:  # unique_values_count == 3
        # Case 3: a, b, and c are all distinct.
        # The maximum value occurs at one vertex and the minimum value at another.
        # The shortest distance is the side length of the triangle.
        distance = side_length
        print("All three vertex values are distinct.")
        print(f"The shortest distance is the triangle's side length: {distance}")

if __name__ == '__main__':
    # You can change these values to test different scenarios.
    # Scenario 1: All values are distinct.
    find_shortest_distance(10, 20, 30)
    print("-" * 20)
    # Scenario 2: Exactly two values are equal.
    find_shortest_distance(15, 15, 5)
    print("-" * 20)
    # Scenario 3: All values are equal.
    find_shortest_distance(7, 7, 7)
