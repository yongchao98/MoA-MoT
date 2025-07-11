import math

def solve_triangle_distance(a, b, c, n):
    """
    Calculates the shortest distance between the nodes with the maximum and minimum values
    in a described equilateral triangle grid.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
        n (int): The number of divisions on each side of the triangle.
    """
    print(f"Solving for a={a}, b={b}, c={c}, n={n}")

    # Case 1: All values are equal
    if a == b and a == c:
        print("All values a, b, c are equal.")
        print("The maximum and minimum values are the same, and all nodes have this value.")
        print("The shortest distance is 0.")
        return 0

    # Case 2: All values are distinct
    if a != b and b != c and a != c:
        print("The values a, b, c are all distinct.")
        print("The maximum and minimum values will be at two different vertices from {A, B, C}.")
        print("The distance between any two vertices of a side length 1 equilateral triangle is 1.")
        print("The shortest distance is 1.")
        return 1

    # Case 3: Exactly two values are equal. The distance may depend on n.
    print("Exactly two of the values a, b, c are equal.")

    def distance_formula(n, i1, j1, i2, j2, print_eq=False):
        di = i1 - i2
        dj = j1 - j2
        # d = (1/n) * sqrt(di^2 + di*dj + dj^2)
        dist_sq_num = di**2 + di * dj + dj**2
        
        if print_eq:
            print(f"The distance is between P(i={i1}, j={j1}) and P(i={i2}, j={j2}).")
            print(f"Using d = (1/n) * sqrt((i1-i2)^2 + (i1-i2)(j1-j2) + (j1-j2)^2)")
            print(f"d = (1/{n}) * sqrt(({di})^2 + ({di})*({dj}) + ({dj})^2)")
            print(f"d = (1/{n}) * sqrt({di**2} + {di*dj} + {dj**2})")
            print(f"d = (1/{n}) * sqrt({dist_sq_num})")

        return math.sqrt(dist_sq_num) / n

    # Find which points hold the max/min
    # There will be one vertex and the opposite side.
    # We need to find the node on the side closest to the vertex.

    # Subcase 3a: Max/min is at A, the other is on side BC.
    # Side BC consists of nodes P(i, n-i) for i in 0..n.
    # We want to minimize distance from A=P(0,0) to P(i, n-i).
    # This minimizes d^2 = i^2 - ni + n^2. Minimum is at i closest to n/2.
    if (b == c and a != b):
        vertex_a_is_max = (a > b)
        print(f"Max/Min value {a} is at vertex A, the other value {b} is on side BC.")
        # Find closest node on BC to A
        i_closest = int(round(n / 2.0))
        p_minmax_on_side = (i_closest, n - i_closest)
        p_vertex = (0, 0)
        
        if vertex_a_is_max:
             max_node, min_node = p_vertex, p_minmax_on_side
        else:
             max_node, min_node = p_minmax_on_side, p_vertex
        dist = distance_formula(n, max_node[0], max_node[1], min_node[0], min_node[1], True)

    # Subcase 3b: Max/min is at B, the other is on side AC.
    # Side AC consists of nodes P(0, j) for j in 0..n.
    # We want to minimize distance from B=P(n,0) to P(0,j).
    # This minimizes d^2 = n^2 - nj + j^2. Minimum is at j closest to n/2.
    elif (a == c and b != a):
        vertex_b_is_max = (b > a)
        print(f"Max/Min value {b} is at vertex B, the other value {a} is on side AC.")
        # Find closest node on AC to B
        j_closest = int(round(n / 2.0))
        p_minmax_on_side = (0, j_closest)
        p_vertex = (n, 0)

        if vertex_b_is_max:
            max_node, min_node = p_vertex, p_minmax_on_side
        else:
            max_node, min_node = p_minmax_on_side, p_vertex
        dist = distance_formula(n, max_node[0], max_node[1], min_node[0], min_node[1], True)

    # Subcase 3c: Max/min is at C, the other is on side AB.
    # Side AB consists of nodes P(i, 0) for i in 0..n.
    # We want to minimize distance from C=P(0,n) to P(i,0).
    # This minimizes d^2 = i^2 - ni + n^2. Minimum is at i closest to n/2.
    elif (a == b and c != a):
        vertex_c_is_max = (c > a)
        print(f"Max/Min value {c} is at vertex C, the other value {a} is on side AB.")
        # Find closest node on AB to C
        i_closest = int(round(n / 2.0))
        p_minmax_on_side = (i_closest, 0)
        p_vertex = (0, n)
        
        if vertex_c_is_max:
            max_node, min_node = p_vertex, p_minmax_on_side
        else:
            max_node, min_node = p_minmax_on_side, p_vertex
        dist = distance_formula(n, max_node[0], max_node[1], min_node[0], min_node[1], True)
        
    print(f"The shortest distance is {dist}")
    return dist


# Example Usage:
# You can change the values of a, b, c, and n here to test different scenarios.
a_val = 5
b_val = 1
c_val = 1
n_val = 3 # Try with an odd number
solve_triangle_distance(a_val, b_val, c_val, n_val)

print("\n" + "="*20 + "\n")

a_val = 5
b_val = 5
c_val = 1
n_val = 4 # Try with an even number
solve_triangle_distance(a_val, b_val, c_val, n_val)