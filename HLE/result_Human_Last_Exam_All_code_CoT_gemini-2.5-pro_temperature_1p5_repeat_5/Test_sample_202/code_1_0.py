import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the node with the maximum value
    and the node with the minimum value in the triangle grid.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    # The value at any node is a weighted average of a, b, and c.
    # The maximum and minimum values in the entire grid are therefore
    # max(a, b, c) and min(a, b, c), respectively.
    # The shortest distance depends on the locations of these extreme values,
    # which is determined by how many of the values a, b, c are unique.

    # We use a set to find the number of unique values.
    unique_values = len(set([a, b, c]))

    distance = 0.0

    if unique_values == 1:
        # Case 1: a = b = c.
        # The maximum and minimum values are the same. All nodes have the same value.
        # The distance is 0.
        distance = 0.0
    elif unique_values == 2:
        # Case 2: Two of the values are equal, one is different.
        # e.g., a = b > c. Max value is on the edge AB, min value is at vertex C.
        # The shortest distance is the altitude of an equilateral triangle with side 1.
        # Altitude = sqrt(side^2 - (side/2)^2) = sqrt(1 - 0.25) = sqrt(0.75) = sqrt(3)/2.
        distance = math.sqrt(3) / 2
    elif unique_values == 3:
        # Case 3: All values a, b, c are distinct.
        # The maximum value is at one vertex and the minimum value is at another.
        # The distance between any two vertices of an equilateral triangle with side 1 is 1.
        distance = 1.0

    print(distance)

# Example Usage:
# You can change these values to test different cases.
# Case 1 (a=b=c): a=5, b=5, c=5 -> Distance: 0
# Case 2 (two equal): a=10, b=10, c=5 -> Distance: sqrt(3)/2
# Case 3 (all distinct): a=10, b=20, c=30 -> Distance: 1
val_a = 10
val_b = 10
val_c = 5

find_shortest_distance(val_a, val_b, val_c)