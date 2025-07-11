import math

def find_shortest_distance(a, b, c):
    """
    Calculates and prints the shortest distance between the node with the maximum
    number and the node with the minimum number in the triangle grid.

    The side length of the equilateral triangle ABC is 1. The numbers a, b, c are
    placed on vertices A, B, C respectively.
    """
    # Identify the number of distinct values among a, b, and c.
    unique_values_count = len(set([a, b, c]))

    print(f"Analyzing for vertex values a={a}, b={b}, c={c}:")

    if unique_values_count == 1:
        # Case 1: a = b = c
        # All nodes have the same value. The maximum and minimum are identical.
        # The shortest distance between the set of max points and min points is 0.
        distance = 0.0
        print("All vertex values are equal.")
        print("The maximum and minimum values are the same everywhere.")
        # The final equation and its numbers
        print("Final Equation: D = 0")
        print(f"Shortest distance: {distance}")

    elif unique_values_count == 2:
        # Case 2: Two values are equal, one is different.
        # One extremum (max or min) is at a single vertex, while the other
        # is along the entire opposite edge. The shortest distance is the
        # altitude of the equilateral triangle with side length 1.
        distance = math.sqrt(3) / 2
        print("Two of the three vertex values are equal.")
        print("The extrema are located at one vertex and the opposite edge.")
        # The final equation and its numbers
        print("Final Equation: D = sqrt(3) / 2")
        print(f"Shortest distance: {distance}")

    else:  # unique_values_count == 3
        # Case 3: a, b, c are all distinct.
        # The maximum value is at one vertex and the minimum is at another.
        # The shortest distance is between these two vertices, i.e., the side length.
        distance = 1.0
        print("All three vertex values are distinct.")
        print("The maximum is at one vertex and the minimum is at another.")
        # The final equation and its numbers
        print("Final Equation: D = 1")
        print(f"Shortest distance: {distance}")

# --- Example Execution ---
# You can change the values of a, b, and c below to test different scenarios.
# The value 'n' (the number of divisions) does not affect the final distance.

print("--- Example 1: All values different ---")
find_shortest_distance(5, 10, 15)
print("\n" + "="*40 + "\n")

print("--- Example 2: Two values equal ---")
find_shortest_distance(10, 10, 5)
print("\n" + "="*40 + "\n")

print("--- Example 3: All values equal ---")
find_shortest_distance(7, 7, 7)