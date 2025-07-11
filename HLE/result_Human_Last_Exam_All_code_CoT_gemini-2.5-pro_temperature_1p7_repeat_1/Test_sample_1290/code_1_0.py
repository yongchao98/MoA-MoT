def solve_maximum_r_vertices():
    """
    Solves for the maximum number of 'r' vertices in ]0, 1[ based on a logical
    deduction from the properties of a simple dessin.

    The argument proceeds in these steps:
    1.  The condition of the dessin being "simple" is analyzed. The key insight
        is that condition (ii) forbids any 'p' or 'q' vertices on the real
        axis from having valency 4. A real vertex with valency 4, being a
        simple extremum, would have two real and two non-real neighbours,
        violating the "all neighbours are real" rule.

    2.  This prohibition has a strong consequence for the poles ('r'-vertices).
        If we assume there are two poles between any two real p/q-vertices,
        the graph of the function phi(x) must have a local extremum between
        those two poles. Since phi is a Belyi function, the value at this
        extremum must be 0 or 1. This would create a real p/q-vertex with
        valency 4, which is a contradiction.

    3.  Therefore, there can be at most one pole ('r'-vertex) between any
        two real 'p' or 'q' vertices.

    4.  The specified form of the function, phi(x) = x^alpha * (1 - x)^beta * F(x),
        implies that x=0 and x=1 are themselves real vertices of the dessin.

    5.  Combining steps 3 and 4, we conclude that there can be at most one
        pole in the interval (0, 1).
    """

    # The maximum number of 'r' vertices is derived from the logical argument.
    max_number_of_r_vertices = 1

    # The final equation is simply setting the maximum number to our derived value.
    print("The final equation for the maximum number of r-vertices (r_max) is:")
    print(f"r_max = {max_number_of_r_vertices}")


solve_maximum_r_vertices()