def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper after a series of folds and cuts.
    """
    print("Let's calculate the total number of edges step by step.")

    # Step 1: Calculate the edges on the outer perimeter.
    # The original 4 corners are cut, adding 4 edges.
    # The original 4 sides each get 2 V-shaped notches. A side with 2 notches is composed of 5 edge segments.
    # For example: ---V---V---
    num_original_sides = 4
    edges_per_notched_side = 5
    side_edges = num_original_sides * edges_per_notched_side
    
    # The 4 corner cuts add 4 more edges connecting the notched sides.
    corner_cut_edges = 4
    
    perimeter_edges = side_edges + corner_cut_edges
    
    print(f"\nFirst, we calculate the edges on the outer perimeter.")
    print(f"Each of the {num_original_sides} sides is notched twice, resulting in {edges_per_notched_side} segments per side.")
    print(f"This gives {num_original_sides} * {edges_per_notched_side} = {side_edges} edges.")
    print(f"The {corner_cut_edges} cuts at the original corners add another {corner_cut_edges} edges.")
    print(f"Total perimeter edges = {side_edges} + {corner_cut_edges} = {perimeter_edges}")

    # Step 2: Calculate the edges from the internal holes.
    # One of the cuts on the folded square creates holes in the interior of the paper.
    # This cut corresponds to the center of the four quadrants of the original paper, creating 4 holes.
    num_holes = 4
    
    # The cut is at a corner formed by two folds. When unfolded, this creates a 4-sided (diamond) hole.
    edges_per_hole = 4
    
    hole_edges = num_holes * edges_per_hole
    
    print(f"\nNext, we calculate the edges from the internal holes.")
    print(f"The cuts create {num_holes} separate holes in the paper.")
    print(f"Each hole has {edges_per_hole} edges.")
    print(f"Total hole edges = {num_holes} * {edges_per_hole} = {hole_edges}")

    # Step 3: Calculate the total number of edges.
    total_edges = perimeter_edges + hole_edges
    
    print(f"\nFinally, we sum the perimeter edges and the hole edges.")
    print(f"Total edges = {perimeter_edges} + {hole_edges} = {total_edges}")

solve_paper_folding_problem()
<<<40>>>