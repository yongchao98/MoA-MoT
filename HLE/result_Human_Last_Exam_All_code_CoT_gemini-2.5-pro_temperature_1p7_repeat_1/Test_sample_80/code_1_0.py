def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard when the
    maximum number of non-attacking bishops are placed on the edge.
    """
    side_length = 8

    # Step 1: Calculate the total number of edge squares.
    # An 8x8 board has 8 squares per side. The four sides have 8*4=32 squares,
    # but the 4 corner squares are counted twice. So, 32 - 4 = 28.
    total_edge_squares = (side_length * 4) - 4

    # Step 2 & 3: Determine the maximum number of non-attacking bishops on the edge.
    # To maximize the number of non-attacking bishops, they must all be placed on
    # squares of the same color. The 28 edge squares are split evenly: 14 white
    # and 14 black. So, the maximum number of bishops is 14.
    max_bishops_placed = total_edge_squares // 2

    # Step 4: Calculate how many edge squares lack bishops.
    # This is the total number of edge squares minus the number occupied by bishops.
    empty_edge_squares = total_edge_squares - max_bishops_placed

    print("Total number of edge squares on an 8x8 board:", total_edge_squares)
    print("Maximum non-attacking bishops that can be placed on the edge:", max_bishops_placed)
    print("Therefore, the number of edge squares that would lack bishops is:")
    print(f"{total_edge_squares} - {max_bishops_placed} = {empty_edge_squares}")

solve_bishop_puzzle()