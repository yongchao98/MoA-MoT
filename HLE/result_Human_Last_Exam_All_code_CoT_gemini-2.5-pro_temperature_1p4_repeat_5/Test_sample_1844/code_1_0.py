def solve_paper_folding_puzzle():
    """
    Calculates the total number of edges on a piece of paper after
    a specific folding, cutting, and unfolding process.
    """
    
    # Part 1: Edges on the outer perimeter
    # The original square has 4 corners. A cut is made on the folded corner
    # that corresponds to the stack of the original 4 corners.
    # This single cut action creates 4 separate notches on the final unfolded shape.
    # Each notch replaces a corner (0 edges) with a new edge (1 edge).
    # The 4 original edges are still present, just shortened.
    # So the total perimeter edges become 4 (original) + 4 (new) = 8.
    perimeter_edges = 8
    
    # Part 2: Edges from the central hole
    # One corner of the folded square corresponds to the exact center of the original paper.
    # A cut here is internal. As we unfold the last two folds, the cut is mirrored
    # twice, creating a single diamond-shaped (4-sided) hole.
    center_hole_edges = 4
    
    # Part 3: Edges from the holes created by the two 'mixed' corners
    # There are two corners on the folded square where one edge is a fold
    # and the other is an original edge of the paper.
    # A cut at one of these 'mixed' corners results in 2 distinct 4-sided holes
    # on the unfolded paper.
    # Since there are 2 such corners, this process creates 2 * 2 = 4 holes in total.
    num_mixed_corners = 2
    holes_per_mixed_cut = 2
    edges_per_hole = 4
    mixed_cuts_edges = num_mixed_corners * holes_per_mixed_cut * edges_per_hole
    
    # Part 4: Sum all the edges
    # The total number of edges is the sum of the perimeter edges and all internal hole edges.
    total_edges = perimeter_edges + center_hole_edges + mixed_cuts_edges
    
    print("This problem can be broken down into the edges from three types of cuts:")
    print(f"1. The outer perimeter of the final shape will have {perimeter_edges} edges.")
    print(f"2. The cut at the center of the paper will create one hole with {center_hole_edges} edges.")
    print(f"3. The other two cuts will create a total of {mixed_cuts_edges} edges from internal holes.")
    print("-" * 20)
    print("The final equation for the total number of edges is:")
    print(f"{perimeter_edges} + {center_hole_edges} + {mixed_cuts_edges} = {total_edges}")

solve_paper_folding_puzzle()
<<<28>>>