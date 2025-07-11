import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the node with the maximum value and
    the node with the minimum value on the described triangular grid.

    The problem states that for any small rhombus in the grid, the sum of the
    numbers on opposite vertices is equal. This property implies that the value
    at any node is a linear interpolation of the values at the main vertices A, B, and C.

    As a result, the maximum and minimum values on the entire grid must occur
    at the vertices A, B, or C themselves. The maximum value will be max(a, b, c)
    and the minimum value will be min(a, b, c).

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    
    Returns:
        float: The shortest distance.
    """
    
    print(f"Analyzing for vertex values a={a}, b={b}, c={c}")
    
    # The maximum and minimum values are determined by the values at the vertices A, B, C.
    max_val = max(a, b, c)
    min_val = min(a, b, c)
    
    print(f"The maximum value on any node is: {max_val}")
    print(f"The minimum value on any node is: {min_val}")

    # Case 1: The max and min values are the same, which means a=b=c.
    # In this case, every node has the same value. We can choose the same node
    # for both max and min, so the distance is 0.
    if max_val == min_val:
        distance = 0.0
        print("All vertex values are equal. The distance is 0.")
    # Case 2: The max and min values are different.
    # This means the nodes with the max and min values must be distinct vertices
    # of the large equilateral triangle (side length 1). The distance
    # between any two distinct vertices is 1.
    else:
        distance = 1.0
        print("The maximum and minimum values occur at distinct vertices of the triangle.")
        print("The distance between any two distinct vertices is 1.")
        
    print(f"Final shortest distance: {distance}\n")
    return distance

if __name__ == '__main__':
    # Example 1: All values are distinct.
    find_shortest_distance(a=10, b=20, c=5)
    
    # Example 2: Two values are the same (but not all three).
    find_shortest_distance(a=15, b=15, c=30)
    
    # Example 3: All values are equal.
    find_shortest_distance(a=7, b=7, c=7)