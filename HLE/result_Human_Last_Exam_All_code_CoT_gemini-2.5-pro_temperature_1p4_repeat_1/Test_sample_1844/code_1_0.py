def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper
    after a series of folds and cuts.
    """

    # Part 1: Calculate the number of edges on the outer boundary.

    # We start with a square, which has 4 edges.
    initial_edges = 4
    print(f"The original square piece of paper has {initial_edges} edges.")
    print("-" * 20)

    # The 4 cuts on the folded square behave differently upon unfolding.
    # Cut 1: At a true corner of the original paper.
    # This single cut clips one of the paper's original 4 corners.
    # Clipping a corner replaces a vertex with an edge, increasing the total edge count by 1.
    corner_clips = 1
    edge_increase_from_corner_clip = 1
    print(f"One of the cuts clips a single corner of the original paper.")
    print(f"This increases the number of outer edges by {edge_increase_from_corner_clip}.")
    
    # Cut 2 & 3: At corners that lie on an original edge of the paper.
    # Because these corners are on a fold, each cut creates two V-shaped notches on the unfolded paper.
    cuts_creating_notches = 2
    notches_per_cut = 2
    # Each V-notch replaces one straight edge segment with two new segments, for a net increase of 1 edge.
    edge_increase_per_notch = 1
    total_notches = cuts_creating_notches * notches_per_cut
    edge_increase_from_notches = total_notches * edge_increase_per_notch
    print(f"\nTwo other cuts create {total_notches} V-shaped notches on the edges of the original paper.")
    print(f"Each notch increases the outer edge count by {edge_increase_per_notch}, adding a total of {edge_increase_from_notches} edges.")

    # Calculate the total number of edges on the outer boundary.
    outer_edges = initial_edges + edge_increase_from_corner_clip + edge_increase_from_notches
    print(f"\nThe final outer boundary will have {initial_edges} + {edge_increase_from_corner_clip} + {edge_increase_from_notches} = {outer_edges} edges.")
    print("-" * 20)

    # Part 2: Calculate the number of edges from internal holes.

    # Cut 4: At the most central corner of the folded paper (a corner formed by two folds).
    # When unfolded, this single cut creates multiple separate holes.
    # The two perpendicular folds cause the cut to be mirrored twice, creating 4 identical holes.
    holes_created_by_central_cut = 4
    # A straight-line cut across a corner formed by two perpendicular folds creates a diamond-shaped (rhombus) hole when unfolded.
    edges_per_hole = 4
    print(f"The fourth cut is at the center and creates {holes_created_by_central_cut} internal holes.")
    print(f"Each hole is diamond-shaped, which has {edges_per_hole} edges.")

    # Calculate the total number of edges for all the holes.
    hole_edges = holes_created_by_central_cut * edges_per_hole
    print(f"The total number of edges from these internal holes is {holes_created_by_central_cut} * {edges_per_hole} = {hole_edges}.")
    print("-" * 20)

    # Part 3: Calculate the total number of edges.
    total_edges = outer_edges + hole_edges
    print(f"To find the total number of edges, we sum the outer and inner edges:")
    print(f"Total Edges = {outer_edges} (outer) + {hole_edges} (holes) = {total_edges}")

solve_paper_folding_problem()