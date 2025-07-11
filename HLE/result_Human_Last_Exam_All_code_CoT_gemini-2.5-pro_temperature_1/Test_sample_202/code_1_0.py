import math

def solve_triangle_distance(n, a, b, c):
    """
    Calculates the shortest distance between the max and min value nodes in the triangle.

    Args:
        n (int): The number of divisions on each side of the triangle.
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    # Using a skewed coordinate system with A at origin (0,0), B at (n,0) and C at (0,n),
    # the value f(i,j) at a node (i,j) can be expressed as a linear function.
    # f(i, j) = k1*i + k2*j + k3
    # We find k1, k2, k3 using the values at vertices A, B, and C.
    k1 = (b - a) / n
    k2 = (c - a) / n
    k3 = a

    print("The value at any node (i, j) is given by the equation:")
    # The prompt requires printing each number in the final equation.
    print(f"f(i, j) = (({b} - {a}) / {n}) * i + (({c} - {a}) / {n}) * j + {a}")
    print(f"f(i, j) = ({k1}) * i + ({k2}) * j + {k3}")
    print("-" * 30)

    # The maximum and minimum of this linear function over the triangular grid must
    # occur at the vertices of the large triangle: A, B, or C.
    values = {'A': a, 'B': b, 'C': c}
    max_val = max(values.values())
    min_val = min(values.values())

    # Identify the vertices corresponding to the max and min values.
    max_points = [point for point, val in values.items() if val == max_val]
    min_points = [point for point, val in values.items() if val == min_val]

    print(f"The maximum value in the grid is: {max_val}")
    print(f"This value is found at vertex/vertices: {', '.join(max_points)}")
    print(f"\nThe minimum value in the grid is: {min_val}")
    print(f"This value is found at vertex/vertices: {', '.join(min_points)}")
    print("-" * 30)

    # Calculate the shortest distance between a max point and a min point.
    if max_val == min_val:
        # This occurs when a = b = c. The max and min points can be the same.
        # The shortest distance is 0.
        distance = 0.0
    else:
        # The max and min values are at different vertices of the triangle ABC.
        # The distance between any two distinct vertices of an equilateral triangle
        # with side length 1 is 1.
        distance = 1.0

    print(f"The shortest distance between a point with the maximum number and a point with the minimum number is: {distance}")
    return distance

if __name__ == '__main__':
    # Example values for n, a, b, c
    n_divs = 10
    val_a = 5.0
    val_b = 20.0
    val_c = 1.0

    solve_triangle_distance(n_divs, val_a, val_b, val_c)
