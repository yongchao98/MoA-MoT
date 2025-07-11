def solve_paper_folding_problem():
    """
    This function calculates the total number of edges on a piece of paper
    after a specific sequence of folding and cutting.
    """
    
    # 1. Calculate the edges from internal holes.
    # A cut at the 'central' corner of the folded square creates 4 holes.
    # Each hole is four-sided.
    num_holes = 4
    edges_per_hole = 4
    internal_edges = num_holes * edges_per_hole

    # 2. Calculate the edges on the outer boundary.
    # The final shape has two sides that are un-notched and two sides with two notches each.
    # An un-notched side with clipped corners has 3 edges.
    # A side with two notches and clipped corners has 9 edges.
    num_unnotched_sides = 2
    edges_per_unnotched_side = 3
    
    num_notched_sides = 2
    edges_per_notched_side = 9

    outer_boundary_edges = (num_unnotched_sides * edges_per_unnotched_side) + (num_notched_sides * edges_per_notched_side)

    # 3. Calculate the total number of edges.
    total_edges = internal_edges + outer_boundary_edges

    # 4. Print the final equation.
    print("The final shape has edges from internal holes and from the outer boundary.")
    print(f"Number of internal edges = {num_holes} holes * {edges_per_hole} edges/hole = {internal_edges}")
    print(f"Number of outer boundary edges = ({num_unnotched_sides} sides * {edges_per_unnotched_side} edges/side) + ({num_notched_sides} sides * {edges_per_notched_side} edges/side) = {outer_boundary_edges}")
    print(f"Total edges = {internal_edges} + {outer_boundary_edges} = {total_edges}")

solve_paper_folding_problem()
<<<40>>>