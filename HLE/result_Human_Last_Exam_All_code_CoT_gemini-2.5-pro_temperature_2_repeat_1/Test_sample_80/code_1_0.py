# The size of the chessboard
n = 8

# Step 1: Calculate the total number of squares on the edge of the board.
# The formula is 4 * n - 4, which accounts for the four corners being counted once.
total_edge_squares = 4 * n - 4

# Step 2: Determine the maximum number of non-attacking bishops that can be placed on the edge squares.
# For an n x n board, the solution is 2 * n - 2.
max_bishops_on_edge = 2 * n - 2

# Step 3: Calculate the number of edge squares that would lack bishops.
# This is the difference between the total edge squares and the number of bishops placed.
empty_squares = total_edge_squares - max_bishops_on_edge

# Print the final result, showing each number in the final equation.
print(f"Total number of edge squares on an {n}x{n} board: {total_edge_squares}")
print(f"Maximum number of non-attacking bishops on the edge: {max_bishops_on_edge}")
print("Number of edge squares lacking bishops is calculated as:")
print(f"{total_edge_squares} - {max_bishops_on_edge} = {empty_squares}")
