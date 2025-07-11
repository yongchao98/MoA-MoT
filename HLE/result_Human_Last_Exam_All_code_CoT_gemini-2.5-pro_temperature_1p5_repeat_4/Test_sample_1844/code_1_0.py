def solve_paper_folding_puzzle():
    """
    Calculates the total number of edges on a piece of paper after a series of folds and cuts.
    """

    # --- Step 1: Analyze the grid ---
    # Folding the paper in half four times creates a grid of 4x4 smaller squares.
    # The grid of squares is 4x4. The grid of vertices (corners of squares) is 5x5.
    
    # --- Step 2: Calculate edges from internal holes ---
    # Holes are formed at the intersections and along the fold lines of the grid.
    
    # At each internal intersection, 4 squares meet. The 4 corner cuts create one 4-sided hole.
    # A 4x4 grid has a 3x3 grid of internal intersections.
    num_internal_intersections = 3 * 3
    edges_per_intersection_hole = 4
    hole_edges_from_intersections = num_internal_intersections * edges_per_intersection_hole
    
    # Along each internal fold line segment, 2 squares meet. Their corner cuts create one 2-sided (V-shaped) hole.
    # A 4x4 grid has 4 horizontal fold lines, each with 3 internal segments, and 4 vertical lines with 3 internal segments.
    num_internal_fold_segments = (4 * 3) + (4 * 3) # 12 vertical + 12 horizontal
    edges_per_fold_line_hole = 2
    hole_edges_from_fold_lines = num_internal_fold_segments * edges_per_fold_line_hole
    
    total_hole_edges = hole_edges_from_intersections + hole_edges_from_fold_lines
    
    # --- Step 3: Calculate edges on the outer boundary ---
    # The outer boundary of the original square is also cut.
    
    # The 4 corners of the original square are clipped, each replaced by 1 edge.
    num_corners = 4
    boundary_edges_from_corner_clips = num_corners * 1
    
    # Each of the 4 sides of the original square has 3 points where a fold line meets it.
    # This creates a 2-edged V-shaped notch.
    num_sides = 4
    notches_per_side = 3
    edges_per_notch = 2
    boundary_edges_from_notches = num_sides * notches_per_side * edges_per_notch
    
    # The original edges are now the straight segments connecting the clips and notches.
    # On each side, there are 4 segments: clip-notch, notch-notch, notch-notch, notch-clip.
    segments_per_side = 4
    boundary_edges_from_segments = num_sides * segments_per_side
    
    total_boundary_edges = boundary_edges_from_corner_clips + boundary_edges_from_notches + boundary_edges_from_segments
    
    # --- Step 4: Calculate total edges and print the result ---
    total_edges = total_hole_edges + total_boundary_edges
    
    print("Let's calculate the total number of edges step-by-step.")
    print("\nFirst, we count the edges of the internal holes:")
    print(f"Holes at 9 internal grid intersections: 9 * 4 = {hole_edges_from_intersections} edges")
    print(f"Holes on 24 internal fold lines: 24 * 2 = {hole_edges_from_fold_lines} edges")
    print(f"Total internal hole edges = {hole_edges_from_intersections} + {hole_edges_from_fold_lines} = {total_hole_edges}")
    
    print("\nNext, we count the edges of the outer boundary:")
    print(f"Edges from 4 clipped corners: 4 * 1 = {boundary_edges_from_corner_clips} edges")
    print(f"Edges from 12 V-shaped notches on the sides: 12 * 2 = {boundary_edges_from_notches} edges")
    print(f"Remaining straight segments of the original boundary: 4 * 4 = {boundary_edges_from_segments} edges")
    print(f"Total outer boundary edges = {boundary_edges_from_corner_clips} + {boundary_edges_from_notches} + {boundary_edges_from_segments} = {total_boundary_edges}")
    
    print("\nFinally, we sum them up:")
    print(f"Total Edges = (Total Hole Edges) + (Total Boundary Edges)")
    print(f"Total Edges = {total_hole_edges} + {total_boundary_edges} = {total_edges}")

solve_paper_folding_puzzle()