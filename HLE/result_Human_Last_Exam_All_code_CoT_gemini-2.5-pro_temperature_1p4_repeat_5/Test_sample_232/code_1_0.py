import sys

def solve():
    """
    Calculates the number of regions an unfolded origami crane paper is divided into.
    The method uses Euler's formula for planar graphs, R = E - V + 1, where R is the
    number of regions, E is the number of edges, and V is the number of vertices.
    """

    # Stage 1: Start with the crease pattern of a standard "preliminary base".
    # This grid consists of a horizontal line, a vertical line, and two diagonals on a square.
    # These lines create a graph with:
    # 9 vertices: 1 at the center, 4 at the corners, and 4 at the midpoints of the edges.
    # 16 edges: 8 segments from the center to the boundary, and 8 segments around the boundary.
    # These divide the square into 8 triangular regions.
    V = 9
    E = 16
    R = E - V + 1
    print(f"Stage 1: Unfolded paper with preliminary base creases.")
    print(f"The creases form a graph with {V} vertices and {E} edges, creating {R} regions.")
    print(f"The equation is: R = E - V + 1")
    print(f"Result: {R} = {E} - {V} + 1 = {E - V + 1}")
    print("-" * 40)

    # Stage 2: Add the folds to turn the preliminary base into a "bird base".
    # This adds 4 "kite" folds and 2 horizontal "limit" folds.
    # - This adds 4 new vertices where the kite folds cross the main diagonals.
    # - This adds 14 new edge segments:
    #   - The 4 main diagonal arms are split by the new vertices (+4 edges).
    #   - The 4 new kite lines are themselves split by those vertices (+8 edges).
    #   - The 2 new horizontal limit lines are added (+2 edges).
    delta_V_stage2 = 4
    delta_E_stage2 = 4 + 8 + 2
    V += delta_V_stage2
    E += delta_E_stage2
    R = E - V + 1
    print(f"Stage 2: After adding bird base creases.")
    print(f"This adds {delta_V_stage2} vertices and {delta_E_stage2} edge segments.")
    print(f"Total Vertices (V) = {V}, Total Edges (E) = {E}")
    print(f"The number of regions is now R = E - V + 1")
    print(f"Result: {R} = {E} - {V} + 1 = {E - V + 1}")
    print("-" * 40)

    # Stage 3: Add the reverse folds for the neck and tail.
    # This adds 4 new crease lines (two for the neck, two for the tail).
    # - This adds 4 new vertices where these new creases terminate on existing kite-fold lines.
    # - This adds 8 new edge segments:
    #   - The 4 new crease lines themselves (+4 edges).
    #   - The 4 existing kite-fold segments that are split by the new vertices (+4 edges).
    delta_V_stage3 = 4
    delta_E_stage3 = 4 + 4
    V += delta_V_stage3
    E += delta_E_stage3
    R = E - V + 1
    print(f"Stage 3: After adding neck and tail reverse folds.")
    print(f"This adds {delta_V_stage3} vertices and {delta_E_stage3} edge segments.")
    print(f"Total Vertices (V) = {V}, Total Edges (E) = {E}")
    print(f"The number of regions is now R = E - V + 1")
    print(f"Result: {R} = {E} - {V} + 1 = {E - V + 1}")
    print("-" * 40)

    # Stage 4: Add the final reverse fold for the head.
    # This adds 1 last crease line on the neck flap.
    # - This adds 2 new vertices where the small head crease begins and ends.
    # - This adds 3 new edge segments:
    #   - The new head-crease line itself (+1 edge).
    #   - The 2 existing neck-crease lines that are split by the new vertices (+2 edges).
    delta_V_stage4 = 2
    delta_E_stage4 = 1 + 2
    V += delta_V_stage4
    E += delta_E_stage4
    R = E - V + 1
    print(f"Stage 4: After adding the final head fold.")
    print(f"This adds {delta_V_stage4} vertices and {delta_E_stage4} edge segments.")
    print(f"Total Vertices (V) = {V}, Total Edges (E) = {E}")
    print(f"The final number of regions is R = E - V + 1")
    print(f"Final Equation: {R} = {E} - {V} + 1")
    
solve()
<<<23>>>