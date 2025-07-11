def find_shortest_distance(a, b, c, n):
    """
    Calculates the shortest distance between the node with the maximum value and
    the node with the minimum value in the described equilateral triangle.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
        n (int): The number of divisions on each side.

    Returns:
        int: The shortest distance, which will be 0 or 1.
    """
    # The reasoning for the solution is as follows:
    # 1. The condition on the rhombuses implies the value at any node (i,j)
    #    is a linear function of its coordinates i and j.
    # 2. A linear function over the triangular domain of nodes will have its
    #    maximum and minimum values at the vertices of that domain.
    # 3. Thus, the maximum value in the grid is max(a, b, c) and the minimum
    #    is min(a, b, c), located at the main vertices A, B, or C.
    # 4. The problem then simplifies to finding the distance between the vertices
    #    holding these max and min values.

    # Check if all vertex values are equal
    if a == b and b == c:
        # If all values are the same, the max and min are the same.
        # The set of points with the max value is the same as the set with the
        # min value, so the shortest distance is 0.
        distance = 0
    else:
        # If the values are not all equal, the max and min are located at
        # different vertices. The side length of the triangle is 1, so the
        # distance between any two distinct vertices is 1.
        distance = 1
    
    return distance

if __name__ == '__main__':
    # You can modify these values to test different scenarios.
    
    # Scenario 1: All vertex values are different
    a1, b1, c1, n1 = 10, 20, 30, 5
    dist1 = find_shortest_distance(a1, b1, c1, n1)
    print(f"Scenario 1: a={a1}, b={b1}, c={c1}, n={n1}")
    print(f"The maximum value is {max(a1,b1,c1)} and the minimum value is {min(a1,b1,c1)}.")
    print(f"The shortest distance is: {dist1}\n")

    # Scenario 2: Two vertex values are the same
    a2, b2, c2, n2 = 15, 5, 15, 10
    dist2 = find_shortest_distance(a2, b2, c2, n2)
    print(f"Scenario 2: a={a2}, b={b2}, c={c2}, n={n2}")
    print(f"The maximum value is {max(a2,b2,c2)} and the minimum value is {min(a2,b2,c2)}.")
    print(f"The shortest distance is: {dist2}\n")

    # Scenario 3: All vertex values are the same
    a3, b3, c3, n3 = 7, 7, 7, 3
    dist3 = find_shortest_distance(a3, b3, c3, n3)
    print(f"Scenario 3: a={a3}, b={b3}, c={c3}, n={n3}")
    print(f"The maximum value is {max(a3,b3,c3)} and the minimum value is {min(a3,b3,c3)}.")
    print(f"The shortest distance is: {dist3}\n")