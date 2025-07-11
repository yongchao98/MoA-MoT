def solve_paper_folding_puzzle():
    """
    This function calculates the total number of edges on a piece of paper
    after a series of folds and cuts.
    """

    # Step 1: Calculate the number of edges from the internal holes.
    # The paper is folded four times (2 top-to-bottom, 2 left-to-right), creating 16 layers.
    # One of the four cuts is at the corner of the folded square that corresponds to the
    # absolute center of the original paper. This cut goes through all 16 layers.
    # When unfolded, this single cut creates 4 separate diamond-shaped holes.
    num_holes = 4
    edges_per_hole = 4
    internal_edges = num_holes * edges_per_hole

    print("--- Calculating Internal Edges (from Holes) ---")
    print(f"The cut at the paper's center creates {num_holes} separate holes when unfolded.")
    print(f"Each hole has {edges_per_hole} edges.")
    print(f"Contribution from internal holes: {num_holes} * {edges_per_hole} = {internal_edges} edges.")
    print("")

    # Step 2: Calculate the number of edges on the outer boundary.
    # The other three cuts affect the paper's original boundary. They clip the 4 original
    # corners and create 2 notches on each of the 4 original sides.
    # We can find the number of outer edges by counting the vertices of the final boundary shape.
    
    # An uncut square has 4 vertices (corners).
    # Clipping the 4 corners replaces each single vertex with 2 new vertices.
    vertices_from_clipping = 4 * 2
    
    # Adding 2 notches to each of the 4 sides adds 1 new vertex per notch.
    num_notches_per_side = 2
    num_sides = 4
    vertices_from_notches = num_notches_per_side * num_sides
    
    # For a single closed shape, the number of edges is equal to the number of vertices.
    outer_edges = vertices_from_clipping + vertices_from_notches

    print("--- Calculating Outer Boundary Edges ---")
    print(f"Clipping the 4 original corners creates {vertices_from_clipping} vertices.")
    print(f"Adding {num_notches_per_side} notches to each of the {num_sides} sides adds {vertices_from_notches} more vertices.")
    print(f"Total vertices on the outer boundary: {vertices_from_clipping} + {vertices_from_notches} = {outer_edges} vertices.")
    print(f"This means the outer boundary has {outer_edges} edges.")
    print("")

    # Step 3: Sum the internal and outer edges for the total.
    total_edges = internal_edges + outer_edges

    print("--- Final Calculation ---")
    print("Total Edges = Internal Edges + Outer Edges")
    print(f"The final equation is: {internal_edges} + {outer_edges} = {total_edges}")

solve_paper_folding_puzzle()
<<<32>>>