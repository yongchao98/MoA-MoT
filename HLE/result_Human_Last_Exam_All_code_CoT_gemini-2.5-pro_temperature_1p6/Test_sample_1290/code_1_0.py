def solve():
    """
    Based on the properties of a simple dessin, the number of poles (r-vertices)
    on the real interval ]0, 1[ can be shown to be at most 1.

    Step 1: A real 'p-vertex' in ]0,1[ must be a local minimum of phi(x), and a real
            'q-vertex' must be a local maximum. This follows from condition (ii) that
            real non-special nodes have real neighbours.

    Step 2: Let's assume there are two poles r_1 and r_2 in ]0,1[. Between them,
            the function phi(x) is continuous. A real rational function must have a
            real critical point between any two real poles.

    Step 3: This critical point 'c' must be a 'p-vertex' or a 'q-vertex'.
            If phi goes from +infinity to +infinity between the poles, the critical
            point is a local minimum. For it to be a q-vertex, phi(c) must be 1.
            This configuration ('q-vertex' as a local minimum) contradicts Step 1.
            If phi goes from -infinity to -infinity, the critical point is a local maximum.
            For it to be a p-vertex, phi(c) must be 0. This configuration
            ('p-vertex' as a local maximum) also contradicts Step 1.

    Step 4: The sign of phi(x) flips at simple poles, so between two consecutive
            poles it would go from +infinity to -infinity (or vice versa), which
            can be monotonic and avoid critical points. However, between any two
            *alternating* poles (e.g. r_1 and r_3), the function goes from +infinity
            to +infinity, forcing the contradiction.

    Step 5: This means we can have at most one pole in the interval ]0, 1[.
            Therefore, the maximum number of r-vertices is 1.
    """
    max_r_vertices = 1
    # The final equation is simply: Maximum number = 1
    # As requested, printing each number in the final equation.
    print(1)

solve()