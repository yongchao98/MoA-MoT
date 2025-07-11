import math

def solve_intersection_problem():
    """
    Calculates the number of intersection points in a hypothetical geometric system.

    In this system, through any point not on a given line, there exist exactly
    three distinct lines parallel to the given line.
    """

    # According to the modified axiom, the number of parallel lines (N) we can draw
    # through a point to another line is 3.
    N = 3
    print(f"Number of parallel lines to a given line through a point (N): {N}\n")

    # We start with a triangle ABC. This gives us 3 vertices and 3 lines (the extended sides).
    # We create 3 groups of new lines based on the modified axiom:
    # 1. Group A: N lines through vertex A, parallel to the opposite side BC.
    # 2. Group B: N lines through vertex B, parallel to the opposite side AC.
    # 3. Group C: N lines through vertex C, parallel to the opposite side AB.
    print("We have 3 groups of new lines, one for each vertex (A, B, C).\n")

    # Now we find the new intersection points.
    # Intersections with the original triangle lines or within the same group of new lines
    # all occur at the original vertices (A, B, or C), which are excluded by the problem.
    # Therefore, new points are only created by the intersection of lines from different groups.

    # Calculate intersections between Group A and Group B.
    # Each of the N lines from Group A intersects each of the N lines from Group B.
    intersections_A_B = N * N
    print(f"Number of intersections between lines from Group A and Group B = {N} * {N} = {intersections_A_B}")

    # Calculate intersections between Group B and Group C.
    intersections_B_C = N * N
    print(f"Number of intersections between lines from Group B and Group C = {N} * {N} = {intersections_B_C}")

    # Calculate intersections between Group C and Group A.
    intersections_C_A = N * N
    print(f"Number of intersections between lines from Group C and Group A = {N} * {N} = {intersections_C_A}\n")

    # The total number of new, distinct points is the sum of these, assuming
    # no three lines (one from each group) are concurrent at a single point.
    total_intersections = intersections_A_B + intersections_B_C + intersections_C_A
    print("Total new intersection points are the sum of these:")
    print(f"{intersections_A_B} + {intersections_B_C} + {intersections_C_A} = {total_intersections}")

solve_intersection_problem()