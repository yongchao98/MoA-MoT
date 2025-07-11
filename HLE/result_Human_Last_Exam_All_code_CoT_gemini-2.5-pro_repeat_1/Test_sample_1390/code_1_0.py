def solve_intersection_problem():
    """
    Calculates the number of new intersection points in a hypothetical geometric system.

    The system's fifth axiom is: "Through any point not on a given line, there exist
    exactly three distinct lines parallel to the given line".

    The calculation is based on a triangle ABC, where we draw all possible parallel
    lines through each vertex to its opposite side.
    """

    # According to the new axiom, this is the number of parallel lines we can draw.
    num_parallels_per_vertex = 3

    # We have three groups of new lines, one for each vertex of the triangle.
    # Group A: Lines through vertex A, parallel to side BC.
    # Group B: Lines through vertex B, parallel to side AC.
    # Group C: Lines through vertex C, parallel to side AB.
    
    print(f"In this system, we can draw {num_parallels_per_vertex} parallel lines through a vertex to the opposite side.")
    print("-" * 50)

    # 1. Calculate intersections between lines from Group A and Group B.
    # Each of the `num_parallels_per_vertex` lines from Group A intersects with each of
    # the `num_parallels_per_vertex` lines from Group B, creating a grid of new points.
    intersections_A_B = num_parallels_per_vertex * num_parallels_per_vertex
    print(f"Number of intersection points from lines through A and lines through B:")
    print(f"Equation: {num_parallels_per_vertex} * {num_parallels_per_vertex} = {intersections_A_B}")
    print("-" * 50)

    # 2. Calculate intersections between lines from Group A and Group C.
    intersections_A_C = num_parallels_per_vertex * num_parallels_per_vertex
    print(f"Number of intersection points from lines through A and lines through C:")
    print(f"Equation: {num_parallels_per_vertex} * {num_parallels_per_vertex} = {intersections_A_C}")
    print("-" * 50)

    # 3. Calculate intersections between lines from Group B and Group C.
    intersections_B_C = num_parallels_per_vertex * num_parallels_per_vertex
    print(f"Number of intersection points from lines through B and lines through C:")
    print(f"Equation: {num_parallels_per_vertex} * {num_parallels_per_vertex} = {intersections_B_C}")
    print("-" * 50)

    # The total number of new points is the sum of these distinct sets of intersections.
    total_intersections = intersections_A_B + intersections_A_C + intersections_B_C

    print("Total number of new intersection points (excluding original vertices A, B, and C):")
    # Here we output each number in the final equation as requested.
    print(f"Final Equation: {intersections_A_B} + {intersections_A_C} + {intersections_B_C} = {total_intersections}")

solve_intersection_problem()