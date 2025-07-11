def solve_paper_folding_puzzle():
    """
    Calculates the total number of edges on a piece of paper after a series of folds and cuts.
    """

    # Step 1: Calculate the number of layers after folding.
    # The paper is folded 4 times.
    num_folds = 4
    layers = 2**num_folds

    # Step 2: Calculate the edges from each type of feature created by the 4 distinct corner cuts.
    # We are calculating the number of edges for each polygon that is formed.

    # Cut 1: The outer boundary.
    # Cutting the 4 corners of the original square turns it into an octagon.
    outer_boundary_edges = 8
    print(f"The outer boundary is transformed into an octagon with {outer_boundary_edges} edges.")

    # Cut 2: The central hole.
    # A cut at the centermost point of the folded paper creates a single large octagonal hole.
    central_hole_edges = 8
    print(f"A single central hole is created with {central_hole_edges} edges.")

    # Cut 3: The axial holes.
    # A cut at another corner creates 4 diamond-shaped holes along the main axes.
    num_axial_holes = 4
    edges_per_axial_hole = 4
    axial_hole_edges = num_axial_holes * edges_per_axial_hole
    print(f"{num_axial_holes} holes are created on the main axes, each with {edges_per_axial_hole} edges, for a total of {axial_hole_edges} edges.")

    # Cut 4: The diagonal holes.
    # The last cut creates 4 smaller V-shaped holes on the diagonals.
    num_diagonal_holes = 4
    edges_per_diagonal_hole = 2
    diagonal_hole_edges = num_diagonal_holes * edges_per_diagonal_hole
    print(f"{num_diagonal_holes} holes are created on the diagonals, each with {edges_per_diagonal_hole} edges, for a total of {diagonal_hole_edges} edges.")

    # Step 3: Sum all the edges.
    total_edges = outer_boundary_edges + central_hole_edges + axial_hole_edges + diagonal_hole_edges
    print("\nCalculating the total number of edges:")
    print(f"Total Edges = {outer_boundary_edges} (outer) + {central_hole_edges} (central) + {axial_hole_edges} (axial) + {diagonal_hole_edges} (diagonal)")
    print(f"Total Edges = {total_edges}")

solve_paper_folding_puzzle()