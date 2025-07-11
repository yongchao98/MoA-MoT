def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper after being
    folded four times, having its corners cut, and then unfolded.
    """

    # Part 1: Calculate the edges of the outer perimeter.
    print("--- Calculating Outer Perimeter Edges ---")

    # Start with a simple square piece of paper.
    edges = 4
    print(f"An initial square has {edges} edges.")

    # The first type of cut is at a corner of the folded square that corresponds
    # to the 4 original corners of the paper. Cutting these turns the 4 vertices
    # of the square into 4 new edges, resulting in an octagon.
    perimeter_edges_after_corner_cuts = 8
    print(f"Cutting the 4 original corners changes the shape to an octagon with {perimeter_edges_after_corner_cuts} edges.")

    # Two other cuts on the folded square correspond to the midpoints of the
    # original paper's edges (2 for horizontal edges, 2 for vertical edges).
    # Unfolding these cuts creates V-shaped notches. Each notch replaces one
    # edge segment with three segments (_ becomes \_/), adding 2 edges to the total count.
    
    # There is a notch on the top, bottom, left, and right sides. That's 4 notches in total.
    num_notches = 4
    edge_increase_per_notch = 2
    
    # We can show the calculation step-by-step for horizontal and vertical notches.
    # The 2 notches on the horizontal edges add 2*2=4 edges.
    perimeter_edges_after_h_notches = perimeter_edges_after_corner_cuts + 2 * edge_increase_per_notch
    print(f"Adding 2 notches on the horizontal edges increases the edge count to {perimeter_edges_after_corner_cuts} + 2 * {edge_increase_per_notch} = {perimeter_edges_after_h_notches}.")
    
    # The 2 notches on the vertical edges add another 2*2=4 edges.
    final_perimeter_edges = perimeter_edges_after_h_notches + 2 * edge_increase_per_notch
    print(f"Adding 2 notches on the vertical edges increases the edge count to {perimeter_edges_after_h_notches} + 2 * {edge_increase_per_notch} = {final_perimeter_edges}.")
    
    print(f"Total outer perimeter edges = {final_perimeter_edges}")
    print("")

    # Part 2: Calculate the edges of the internal holes.
    print("--- Calculating Internal Hole Edges ---")

    # The last cut is at the corner of the folded square that corresponds to the
    # center of the original paper. This corner is formed by two perpendicular creases.
    # Unfolding a cut at such a corner creates a hole.
    # The process of folding twice (T->B, L->R) and then cutting the central corner
    # creates a 4-sided hole (a rhombus).
    # Since this folding process was performed twice, the unfolding process reveals
    # that the rhombus hole itself is duplicated for each quadrant of the original paper.
    num_holes = 4
    edges_per_hole = 4
    total_internal_edges = num_holes * edges_per_hole
    print(f"The central cut creates {num_holes} separate holes.")
    print(f"Each hole is a rhombus, which has {edges_per_hole} edges.")
    print(f"Total internal edges = {num_holes} * {edges_per_hole} = {total_internal_edges}.")
    print("")

    # Part 3: Calculate the total number of edges.
    print("--- Calculating Total Edges ---")
    total_edges = final_perimeter_edges + total_internal_edges
    print("Total edges = Perimeter Edges + Internal Edges")
    print(f"Total edges = {final_perimeter_edges} + {total_internal_edges} = {total_edges}")

solve_paper_folding_problem()