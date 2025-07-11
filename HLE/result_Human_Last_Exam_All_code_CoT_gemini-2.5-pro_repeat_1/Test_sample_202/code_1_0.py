import math

def find_shortest_distance_in_triangle(a, b, c):
    """
    Calculates the shortest distance between the node with the maximum value and
    the node with the minimum value in the described triangle grid.

    The problem describes an equilateral triangle ABC with side length 1.
    The values at the vertices A, B, C are a, b, c.
    The rules for assigning values to the internal nodes result in the value
    at any node being a linear interpolation of a, b, and c.
    
    As a result, the maximum and minimum values in the entire grid must occur
    at the vertices A, B, or C.

    The shortest distance depends on how many of a, b, c are equal.
    - If a, b, c are all distinct, the max is at one vertex and the min is at another.
      The distance is the side length of the triangle, which is 1.
    - If two of a, b, c are equal (e.g., a=b > c), the extremum is an entire side
      and the other is the opposite vertex. The distance is the triangle's altitude.
      For a side length of 1, the altitude is sqrt(3)/2.
    - If a, b, c are all equal, all points have the same value, so the distance is 0.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.

    Returns:
        float: The shortest distance between the max and min valued points.
    """
    
    # Sort the values to easily identify the cases
    values = sorted([a, b, c])
    v_min = values[0]
    v_mid = values[1]
    v_max = values[2]

    # Case 1: All three values are equal.
    # The max and min points are the same set (the whole triangle).
    if v_min == v_max:
        return 0.0

    # Case 2: Exactly two values are equal.
    # This happens if the min two are equal or the max two are equal.
    # The distance is the altitude of the triangle.
    elif (v_min == v_mid) or (v_mid == v_max):
        return math.sqrt(3) / 2

    # Case 3: All three values are distinct.
    # The distance is between two vertices of the triangle.
    else:  # v_min < v_mid < v_max
        return 1.0

# --- Example Usage ---
# You can change these values to test different scenarios.
# The value of `n` is not needed for the calculation.

# Example 1: a, b, c are all different
val_a1, val_b1, val_c1 = 10, 5, 2

# Example 2: Two values are equal (one max, two mins)
val_a2, val_b2, val_c2 = 10, 3, 3

# Example 3: All values are equal
val_a3, val_b3, val_c3 = 7, 7, 7

# Calculate and print the distance for the first example
# To solve your specific problem, replace the values here.
distance = find_shortest_distance_in_triangle(val_a1, val_b1, val_c1)
print(distance)