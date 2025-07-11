def find_shortest_distance():
    """
    Calculates the shortest distance between the node with the maximum number
    and the node with the minimum number in the described triangle grid.
    """
    # --- User-definable parameters ---
    # The values on the vertices A, B, and C of the equilateral triangle.
    a = 10.0
    b = 20.0
    c = 10.0

    # The number of divisions on each side.
    # Note: This value does not affect the final distance.
    n = 10
    # ------------------------------------

    # As derived from the problem's conditions, the value at any node
    # is a linear interpolation of the values at the vertices A, B, C.
    # This means the maximum and minimum values for the entire grid must
    # occur at the vertices A, B, or C.
    max_val = max(a, b, c)
    min_val = min(a, b, c)

    # The equation for the value V at a node with barycentric coordinates (i/n, j/n, k/n)
    # where i+j+k=n is: V = a*(i/n) + b*(j/n) + c*(k/n)
    # The numbers in this equation are:
    print(f"Given values at vertices: a = {a}, b = {b}, c = {c}")
    print(f"Number of divisions per side: n = {n}")
    print("-" * 30)

    # Determine the distance based on whether the max and min values are the same.
    if max_val == min_val:
        # If max equals min, then a = b = c. All nodes have the same value.
        # The set of maximum points and minimum points are identical.
        # The shortest distance between them is 0.
        distance = 0
    else:
        # If max is different from min, the maximum value is located at one
        # vertex and the minimum value at another. The distance between any
        # two vertices of an equilateral triangle with side length 1 is 1.
        distance = 1

    print(f"The maximum value on any node is: {max_val}")
    print(f"The minimum value on any node is: {min_val}")
    print(f"The shortest distance between a max point and a min point is: {distance}")


# Execute the function to find and print the result.
find_shortest_distance()