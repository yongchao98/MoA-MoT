import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the nodes with the maximum and
    minimum values on the triangle grid.

    The problem describes an equilateral triangle with side length 1.
    The value at any node is a linear interpolation of the values a, b, c at
    the main vertices A, B, C. This means the maximum and minimum values
    will always be located on the boundaries of the triangle (vertices or edges).

    Args:
        a (float): The number placed on vertex A.
        b (float): The number placed on vertex B.
        c (float): The number placed on vertex C.
    """
    values = [a, b, c]
    max_val = max(values)
    min_val = min(values)

    print(f"Given the numbers a={a}, b={b}, c={c} at the vertices.")
    print(f"The maximum value on the grid is {max_val}.")
    print(f"The minimum value on the grid is {min_val}.")
    print("-" * 30)

    # Case 1: All values are the same.
    if max_val == min_val:
        print("All three numbers are equal.")
        print("The maximum and minimum values are the same, so the distance is 0.")
        distance = 0.0
        print(f"\nFinal Equation: Distance = {distance}")
        return

    count_max = values.count(max_val)
    count_min = values.count(min_val)

    # Case 2: Two values are equal (and are the max or min).
    # This means one extremum is at a vertex, and the other is along the opposite edge.
    # e.g., a=b > c. Max is on edge AB, min is at vertex C.
    if count_max == 2 or count_min == 2:
        print("Two of the three numbers are equal, forming an extremum.")
        print("One extremum lies on an edge, and the other is at the opposite vertex.")
        print("The shortest distance is the altitude of the equilateral triangle.")
        side = 1.0
        altitude = math.sqrt(side**2 - (side / 2)**2)
        print("\nFinal Equation: Altitude = sqrt(side^2 - (side/2)^2)")
        print(f"Calculation: sqrt({side}^2 - ({side/2})^2) = sqrt({side**2 - (side/2)**2}) = {altitude}")
    # Case 3: All values are distinct.
    # The max is at one vertex, the min is at another.
    else:
        print("The three numbers are distinct.")
        print("The maximum is at one vertex and the minimum is at another.")
        print("The shortest distance is the side length of the triangle.")
        side = 1.0
        print(f"\nFinal Equation: Distance = Side Length = {side}")

# --- You can change the values of a, b, and c below to test different scenarios ---

# Example 1: Distinct values (a=10, b=5, c=1)
# Expected distance: 1
print("--- Scenario 1: Distinct Values ---")
find_shortest_distance(10, 5, 1)
print("\n" + "="*40 + "\n")

# Example 2: Two values are max (a=10, b=10, c=1)
# Expected distance: sqrt(3)/2
print("--- Scenario 2: Two Max Values ---")
find_shortest_distance(10, 10, 1)
print("\n" + "="*40 + "\n")

# Example 3: Two values are min (a=10, b=1, c=1)
# Expected distance: sqrt(3)/2
print("--- Scenario 3: Two Min Values ---")
find_shortest_distance(10, 1, 1)
print("\n" + "="*40 + "\n")

# Example 4: All values equal (a=5, b=5, c=5)
# Expected distance: 0
print("--- Scenario 4: All Values Equal ---")
find_shortest_distance(5, 5, 5)
