# The size of a standard chessboard
board_size = 8

# Step 1: Calculate the total number of edge squares.
# On an N x N board, there are N squares on the top row, N on the bottom row,
# and (N - 2) on both the left and right columns (as corners are already counted).
# Total edge squares = N + N + (N-2) + (N-2) = 4*N - 4.
total_edge_squares = 4 * (board_size - 1)

print(f"First, we calculate the total number of edge squares on a {board_size}x{board_size} board.")
print(f"There are {board_size} squares on each of the 4 sides, but the 4 corner squares are counted twice.")
print(f"So, the total number of edge squares is 4 * {board_size} - 4 = {total_edge_squares}.")
print("-" * 20)

# Step 2: Determine the maximum number of non-attacking bishops on the edge.
# The maximum number of non-attacking bishops on an N x N board is 2*N - 2.
# For an 8x8 board, this is 2 * 8 - 2 = 14.
# This number can be achieved by placing bishops only on edge squares. For example,
# placing bishops on all 8 squares of the first row (a1 to h1) and all 6 squares
# of the eighth row excluding the corners (b8 to g8).
# This gives a total of 8 + 6 = 14 non-attacking bishops.
max_bishops = 2 * board_size - 2

print("Next, we determine the maximum number of bishops that can be placed on these edge squares without any attacking each other.")
print(f"The maximum for a {board_size}x{board_size} board is 2 * {board_size} - 2 = {max_bishops}.")
print("-" * 20)

# Step 3: Calculate the number of edge squares lacking bishops.
empty_edge_squares = total_edge_squares - max_bishops

print("Finally, we subtract the number of bishops from the total number of edge squares to find the number of empty squares.")
# The final equation as requested
print(f"The final calculation is: {total_edge_squares} - {max_bishops} = {empty_edge_squares}")

<<<14>>>