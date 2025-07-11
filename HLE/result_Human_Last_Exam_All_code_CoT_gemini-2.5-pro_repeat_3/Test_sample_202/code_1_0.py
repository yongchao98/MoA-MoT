import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the node with the maximum value
    and the node with the minimum value in the described triangle grid.

    The problem statement implies that the value at any node is a linear
    interpolation of the values at the main vertices A, B, and C.
    The value at a node is given by V = u*a + v*b + w*c, where (u,v,w) are
    the barycentric coordinates of the node.

    The maximum and minimum of this linear function over the triangle are
    always max(a, b, c) and min(a, b, c), respectively. The location of
    these extremal values depends on how many vertices share them.

    Args:
        a: The number on vertex A.
        b: The number on vertex B.
        c: The number on vertex C.
    """
    print(f"Given values at vertices A, B, C: a={a}, b={b}, c={c}\n")

    # We count the number of distinct values among a, b, and c to determine the case.
    unique_values = set([a, b, c])
    num_unique = len(unique_values)

    if num_unique == 3:
        # Case 1: a, b, and c are all different.
        # The maximum value occurs at exactly one vertex.
        # The minimum value occurs at exactly one other vertex.
        # The shortest distance is the distance between these two vertices,
        # which is the side length of the equilateral triangle.
        distance = 1.0
        print("The three vertex values are distinct.")
        print("The maximum is at one vertex and the minimum is at another.")
        print("The shortest distance is the side length of the triangle.")
        print("Final Answer: The distance is 1.")

    elif num_unique == 2:
        # Case 2: Exactly two of the values are equal.
        # e.g., if a = b > c, the max value is on edge AB and the min is at vertex C.
        # The shortest distance is from a vertex to the opposite side, which is the
        # altitude of the equilateral triangle.
        distance = math.sqrt(3) / 2
        print("Exactly two of the three vertex values are equal.")
        print("One extremal value (max or min) lies on an edge, while the other is at the opposite vertex.")
        print("The shortest distance is the altitude of the triangle.")
        # The final equation is Distance = sqrt(3) / 2
        print(f"Final Answer: The distance is sqrt(3) / 2 = {distance}")

    else: # num_unique == 1
        # Case 3: All three values are equal, a = b = c.
        # The maximum and minimum values are the same and occur everywhere.
        # The shortest distance is 0.
        distance = 0.0
        print("All three vertex values are equal.")
        print("The maximum and minimum values are the same, so the sets of points are identical.")
        print("The shortest distance is 0.")
        print("Final Answer: The distance is 0.")

# --- Example Usage ---
# You can change these values to test different scenarios.
# Scenario 1: Three distinct values (e.g., 5, 10, 15) -> distance = 1
# Scenario 2: Two distinct values (e.g., 5, 5, 10) -> distance = sqrt(3)/2
# Scenario 3: One distinct value (e.g., 8, 8, 8) -> distance = 0
val_a = 5
val_b = 5
val_c = 10

find_shortest_distance(val_a, val_b, val_c)