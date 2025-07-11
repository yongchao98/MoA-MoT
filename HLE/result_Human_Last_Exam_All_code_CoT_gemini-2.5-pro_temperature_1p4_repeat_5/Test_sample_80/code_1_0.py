def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing the maximum
    number of non-attacking bishops on the edge.
    """
    
    # A standard chessboard has 8 squares per side.
    side_length = 8
    
    # 1. Calculate the total number of edge squares.
    # Formula: 4 sides * length_of_side - 4 corners (to avoid double-counting).
    total_edge_squares = (side_length * 4) - 4
    
    # 2. Determine the maximum number of non-attacking bishops on the edge.
    # A good strategy is to place bishops on the two side files ('a' and 'h').
    # This gives 8 squares on file 'a' and 8 on file 'h'.
    squares_on_side_files = side_length * 2
    
    # However, bishops on a1 and h8 attack each other, and bishops on a8 and h1 attack each other.
    # This creates 2 pairs of mutually attacking bishops.
    attacking_pairs = 2
    
    # To have a valid placement, we must remove one bishop from each attacking pair.
    # So, the maximum number of bishops we can place is:
    max_bishops = squares_on_side_files - attacking_pairs
    
    # 3. Calculate the number of edge squares without bishops.
    empty_squares = total_edge_squares - max_bishops
    
    # 4. Print the final equation showing each number.
    print(f"Total edge squares: {total_edge_squares}")
    print(f"Maximum non-attacking bishops on the edge: {max_bishops}")
    print("The number of edge squares that would lack bishops is:")
    print(f"{total_edge_squares} - {max_bishops} = {empty_squares}")

solve_bishop_puzzle()