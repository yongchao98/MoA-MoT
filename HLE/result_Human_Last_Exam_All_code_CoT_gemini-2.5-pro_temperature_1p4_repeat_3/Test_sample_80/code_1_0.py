# Step 1: Define the size of the chessboard.
board_size = 8

# Step 2: Calculate the total number of edge squares on the board.
# Formula: (squares per side * 4 sides) - 4 corners (to avoid double-counting).
total_edge_squares = (board_size * 4) - 4

# Step 3: Calculate the maximum number of non-attacking bishops that can be placed on the edge.
# On an N x N board (where N is even), the edge squares of each color can be formed into pairs
# where the squares in each pair attack each other. Therefore, you can only place one bishop
# per pair. The maximum number of bishops is half the total number of edge squares.
max_bishops = total_edge_squares // 2

# Step 4: Calculate the number of edge squares that will remain empty.
empty_squares = total_edge_squares - max_bishops

# Step 5: Print the final equation and the answer.
print("This problem can be solved by first finding the total number of edge squares,")
print("then finding the maximum number of non-attacking bishops that can be placed on them.")
print("The number of empty squares is the difference between these two values.")
print("\nEquation:")
print(f"Total Edge Squares ({total_edge_squares}) - Max Bishops ({max_bishops}) = Empty Squares ({empty_squares})")
print(f"\nSo, {empty_squares} edge squares would lack bishops.")
