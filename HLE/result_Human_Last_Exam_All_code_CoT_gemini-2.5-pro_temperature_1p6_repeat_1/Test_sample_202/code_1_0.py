import math

def solve_triangle_distance(a, b, c):
    """
    Calculates the shortest distance between the node with the maximum value
    and the node with the minimum value in the described triangle problem.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    
    # The side length of the main equilateral triangle is 1.
    side_length = 1.0
    
    # Store the vertex values in a list
    values = [a, b, c]
    
    # Find the maximum and minimum values
    v_max = max(values)
    v_min = min(values)
    
    # Print the problem setup
    print(f"Vertex values are a={a}, b={b}, c={c}")
    print(f"Maximum value on the grid: {v_max}")
    print(f"Minimum value on the grid: {v_min}")
    
    # Case 1: All values are the same. Max and min are equal.
    if v_max == v_min:
        distance = 0.0
        print("All vertices have the same value.")
        print("The distance is 0.")
        return

    # Count how many vertices have the max value and how many have the min value.
    count_max = values.count(v_max)
    count_min = values.count(v_min)
    
    # Case 2: Max and min values are on single, distinct vertices.
    # This happens when all three values a, b, c are different.
    if count_max == 1 and count_min == 1:
        distance = side_length
        print("Max and Min values are at distinct vertices.")
        # The final equation is simply the distance being the side length.
        print(f"The shortest distance is {int(distance)}.")

    # Case 3: One of the extreme values is shared by two vertices.
    # The set of points for one extreme value is an edge, and the other extreme is at the opposite vertex.
    else:
        # The distance is the altitude of the equilateral triangle.
        # h = sqrt(side^2 - (side/2)^2) = sqrt(1 - 1/4) = sqrt(3/4) = sqrt(3)/2
        distance = math.sqrt(3) / 2
        print("One extreme value is on an edge, the other is at the opposite vertex.")
        # Final equation for the distance
        print(f"The shortest distance is sqrt(3)/2 = {distance}")

# Example usage with values where the distance is the altitude of the triangle.
# For example, a=5, b=5, c=2. The max value 5 is on edge AB. The min value 2 is at C.
# The shortest distance is from C to edge AB.
solve_triangle_distance(5, 5, 2)
