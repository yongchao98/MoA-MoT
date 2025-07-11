import math

def find_shortest_distance(a, b, c, n):
    """
    Calculates the shortest distance between the nodes with the maximum and
    minimum values on the described equilateral triangle.

    The problem states that in any small rhombus, the sum of opposite vertices
    is equal. This implies that the value at any node is a linear function
    of its position. Therefore, the maximum and minimum values for the entire
    set of nodes must occur at the main vertices A, B, or C.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
        n (int): The number of divisions on each side of the triangle.
    """

    print(f"Analyzing for values: a = {a} (at vertex A), b = {b} (at vertex B), c = {c} (at vertex C).")

    # The maximum value in the triangle will be max(a, b, c)
    max_val = max(a, b, c)
    # The minimum value in the triangle will be min(a, b, c)
    min_val = min(a, b, c)
    
    print(f"The maximum value on any node is {max_val}.")
    print(f"The minimum value on any node is {min_val}.")

    # Case 1: All three vertices have the same value.
    if a == b and b == c:
        print("All vertices have the same value. The max and min values are identical.")
        distance = 0.0
        print(f"Therefore, the shortest distance is: {distance}")
        return distance

    # Count how many vertices share the max and min values
    values = [a, b, c]
    max_count = values.count(max_val)
    min_count = values.count(min_val)

    # Case 2: Two vertices share an extreme value.
    # The extremum is on an edge, and the other is at the opposite vertex.
    # The distance is the altitude of the triangle.
    if max_count == 2 or min_count == 2:
        if max_count == 2:
             print(f"The maximum value is shared by two vertices.")
             print(f"The minimum value is at a single vertex.")
        else: # min_count == 2
             print(f"The minimum value is shared by two vertices.")
             print(f"The maximum value is at a single vertex.")
        
        distance = math.sqrt(3) / 2
        print("The shortest distance is between a vertex and the opposite side.")
        print(f"This distance is the altitude of the triangle: sqrt(3)/2 ≈ {distance:.4f}")
        return distance

    # Case 3: All three values are distinct.
    # The max and min values are at two different vertices.
    else:
        print("The maximum and minimum values are at two distinct vertices.")
        distance = 1.0
        print(f"The shortest distance is the side length between these two vertices, which is: {distance}")
        return distance

# --- Example Execution ---
# You can change these values to test different scenarios.
# Let's use a case where the values are distinct.
A_val = 5
B_val = 12
C_val = 10
N_divisions = 100

final_distance = find_shortest_distance(A_val, B_val, C_val, N_divisions)
# The prompt requests the final answer in a specific format,
# but since the answer depends on the inputs a, b, c, the script above is the general solution.
# For the specific example executed here (a=5, b=12, c=10), the answer is 1.0.
# The code correctly outputs this result.