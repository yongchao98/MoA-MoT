# Define the dimensions of the chessboard.
rows = 7
cols = 8

# The strategy is to remove all pieces from a central row and a central column.
# Let's choose row 4 and column 4.
num_pieces_in_row_4 = cols
num_pieces_in_col_4 = rows

# The piece at the intersection of the chosen row and column is counted in both removals.
# We need to subtract one to avoid double-counting.
intersection = 1

# Calculate the total number of pieces to be removed.
min_removed = num_pieces_in_row_4 + num_pieces_in_col_4 - intersection

# Print the explanation and the result.
print("A 7x8 rectangular chessboard has dimensions {}x{}.".format(rows, cols))
print("To break all possible lines of 5 or more, we can remove all pieces from a central row and a central column.")
print("Let's choose to remove all pieces from one row and one column.")
print("Number of pieces in a row = {}.".format(cols))
print("Number of pieces in a column = {}.".format(rows))
print("The piece at the intersection is counted twice, so we subtract 1.")
print("The total minimum number of pieces to be removed is:")
print("{} + {} - {} = {}".format(cols, rows, intersection, min_removed))