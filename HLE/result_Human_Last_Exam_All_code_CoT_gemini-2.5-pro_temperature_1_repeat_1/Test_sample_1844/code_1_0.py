def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper after a series of folds and cuts.
    """
    # The paper is folded to create a 4x4 grid of squares.
    n = 4
    print(f"The folding process divides the paper into a {n}x{n} grid of squares.")

    # 1. Calculate horizontal edges
    # A grid of n*n squares has (n+1) horizontal lines.
    # Each line is composed of n segments, which become edges.
    num_horizontal_lines = n + 1
    horizontal_edges_per_line = n
    total_horizontal_edges = num_horizontal_lines * horizontal_edges_per_line
    print(f"There are {num_horizontal_lines} horizontal grid lines, each with {horizontal_edges_per_line} segments.")
    print(f"Total horizontal edges = {num_horizontal_lines} * {horizontal_edges_per_line} = {total_horizontal_edges}")

    # 2. Calculate vertical edges
    # Similarly, there are (n+1) vertical lines, each with n segments.
    num_vertical_lines = n + 1
    vertical_edges_per_line = n
    total_vertical_edges = num_vertical_lines * vertical_edges_per_line
    print(f"There are {num_vertical_lines} vertical grid lines, each with {vertical_edges_per_line} segments.")
    print(f"Total vertical edges = {num_vertical_lines} * {vertical_edges_per_line} = {total_vertical_edges}")

    # 3. Calculate diagonal edges from the cuts
    # Cutting the 4 corners of the folded stack is equivalent to cutting all 4 corners of each of the n*n squares.
    # Each cut creates one new diagonal edge.
    num_squares = n * n
    cuts_per_square = 4
    total_diagonal_edges = num_squares * cuts_per_square
    print(f"There are {num_squares} squares in the grid, and each has its {cuts_per_square} corners cut.")
    print(f"Total diagonal edges = {num_squares} * {cuts_per_square} = {total_diagonal_edges}")

    # 4. Calculate total edges
    total_edges = total_horizontal_edges + total_vertical_edges + total_diagonal_edges
    print("\nCalculating the total number of edges:")
    print(f"Total Edges = (Horizontal Edges) + (Vertical Edges) + (Diagonal Edges)")
    print(f"Total Edges = {total_horizontal_edges} + {total_vertical_edges} + {total_diagonal_edges}")
    print(f"The total number of edges on the unfolded shape is: {total_edges}")

solve_paper_folding_problem()
<<<104>>>