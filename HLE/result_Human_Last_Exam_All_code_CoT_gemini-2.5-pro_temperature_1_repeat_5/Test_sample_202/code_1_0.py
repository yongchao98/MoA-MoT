import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the maximum and minimum valued nodes
    on the described triangular grid.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    print(f"--- Analyzing for vertex values a={a}, b={b}, c={c} ---")
    vals = [a, b, c]
    v_max = max(vals)
    v_min = min(vals)

    # Case 1: All vertex values are equal
    if v_max == v_min:
        print("All vertex values are equal.")
        print("The maximum and minimum points are the same, covering the entire triangle.")
        distance = 0
        print(f"The shortest distance is: {distance}")
        return

    count_max = vals.count(v_max)
    count_min = vals.count(v_min)

    print(f"Maximum value: {v_max}, located at {count_max} vertex/vertices.")
    print(f"Minimum value: {v_min}, located at {count_min} vertex/vertices.")

    # Case 2: Unique maximum and unique minimum
    if count_max == 1 and count_min == 1:
        # This means the max is at one vertex and the min is at another.
        # e.g., a > b > c
        # The shortest distance is between two vertices of the triangle.
        distance = 1.0
        print("The max and min values each occur at a single, distinct vertex.")
        print(f"The shortest distance is the side length of the triangle: {distance}")

    # Case 3: One extremum is on a side, the other is at the opposite vertex
    else:
        # This covers when one value is unique and the other two are equal.
        # e.g., a = b > c OR a > b = c
        # The shortest distance is the height of the equilateral triangle.
        side_length = 1
        height = math.sqrt(side_length**2 - (side_length / 2)**2)
        print("One of the extreme values occurs on a side (2 vertices),")
        print("and the other occurs at the opposite vertex.")
        print("The shortest distance is the height of the triangle.")
        # Outputting the numbers in the final equation for clarity
        print(f"h = sqrt({side_length}^2 - ({side_length / 2})^2) = sqrt(3)/2")
        distance = height
        print(f"The shortest distance is: {distance}")

if __name__ == '__main__':
    # You can change these values to test different scenarios
    example_a, example_b, example_c = 10, 10, 3

    find_shortest_distance(example_a, example_b, example_c)

    print("\n--- Other Examples ---")
    find_shortest_distance(10, 5, 1) # Should be 1
    find_shortest_distance(10, 10, 10) # Should be 0
