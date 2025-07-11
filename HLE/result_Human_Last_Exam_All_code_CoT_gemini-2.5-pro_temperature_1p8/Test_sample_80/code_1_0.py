def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on them.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An n x n board has 4 sides of length n, but the 4 corners are counted twice.
    # So, the number of edge squares is (n * 4) - 4.
    total_edge_squares = (board_size * 4) - 4
    
    # Step 2 & 3: Determine the maximum number of bishops for each color.
    # Bishops on light squares don't attack bishops on dark squares.
    # The problem is symmetric for both colors.
    
    # Let's analyze one color, e.g., dark squares.
    # To place the maximum number of non-attacking bishops, we can fill the
    # squares on the outermost files, 'a' and 'h'.
    # Dark edge squares on file 'a': a1, a3, a5, a7 (4 squares)
    # Dark edge squares on file 'h': h2, h4, h6, h8 (4 squares)
    # Total potential spots = 4 + 4 = 8
    
    # However, bishops on a1 and h8 attack each other along the main diagonal.
    # We can only place one bishop on this diagonal. So, we must remove one.
    # Maximum bishops for one color = 8 - 1 = 7.
    max_bishops_per_color = (board_size // 2) + (board_size // 2) - 1
    
    # The number is the same for light squares due to symmetry.
    # (a8 and h1 attack each other for the light squares).
    max_light_bishops = max_bishops_per_color
    max_dark_bishops = max_bishops_per_color
    
    # Step 4: Calculate the total number of bishops that can be placed.
    total_bishops = max_light_bishops + max_dark_bishops
    
    # Step 5: Calculate how many edge squares are left empty.
    empty_edge_squares = total_edge_squares - total_bishops

    # Final Output
    print("Chessboard Bishops on the Edge Puzzle:")
    print("-" * 40)
    print(f"1. Total edge squares on an {board_size}x{board_size} board: {total_edge_squares}")
    print(f"2. Max non-attacking bishops on light edge squares: {max_light_bishops}")
    print(f"3. Max non-attacking bishops on dark edge squares: {max_dark_bishops}")
    print(f"4. Total bishops that can be placed on the edge: {max_light_bishops} + {max_dark_bishops} = {total_bishops}")
    print("-" * 40)
    print("Question: How many edge squares would lack bishops?")
    print(f"Calculation: {total_edge_squares} (total edge squares) - {total_bishops} (placed bishops)")
    print(f"Answer: {empty_edge_squares}")


solve_chess_puzzle()
<<<14>>>