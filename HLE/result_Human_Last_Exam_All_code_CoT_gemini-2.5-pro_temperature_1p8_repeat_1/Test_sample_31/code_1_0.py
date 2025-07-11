# The last layer of a Rubik's Cube has 8 pieces that are permuted during the PLL step.
num_corners = 4
num_edges = 4

# To fully identify the permutation, we must know the position of each of these pieces.
# We can determine the position of a piece by looking at at least one of its non-top-facing stickers.
# Therefore, the minimum number of stickers we need to see is equal to the number of pieces.
total_pieces_to_identify = num_corners + num_edges

# Final calculation
print(f"Number of corner pieces: {num_corners}")
print(f"Number of edge pieces: {num_edges}")
print(f"To uniquely identify the location of every piece in the last layer,")
print(f"one must observe a sticker on each of the {total_pieces_to_identify} pieces.")
print(f"Therefore, the minimum number of non-top-facing stickers that must be seen is {total_pieces_to_identify}.")
