# The size of the chessboard.
ROWS = 7
COLS = 8

# Step 1: Remove an entire row to break all vertical and diagonal lines.
# We choose a central row (e.g., row 4) for this.
# The number of pieces removed is equal to the number of columns.
removals_for_main_row = COLS
print(f"Step 1: To break all vertical and diagonal lines, we remove a full row.")
print(f"This requires removing {removals_for_main_row} pieces.")

# Step 2: Break the remaining horizontal lines.
# After removing one row, we are left with ROWS - 1 rows.
# Each of these rows is a line of 8 pieces that must be broken.
# To break a line of 8, we need to remove at least 1 piece.
remaining_rows = ROWS - 1
removals_for_other_rows = remaining_rows * 1  # 1 removal per remaining row.
print(f"\nStep 2: There are {remaining_rows} rows remaining, each forming a line of 8.")
print(f"To break these horizontal lines, we must remove at least 1 piece from each.")
print(f"This requires an additional {removals_for_other_rows} removals.")

# Step 3: Calculate the total minimum removals.
total_removals = removals_for_main_row + removals_for_other_rows
print(f"\nStep 3: The minimum total number of pieces to remove is the sum of these two steps.")
print("The final equation is:")
print(f"{removals_for_main_row} + {removals_for_other_rows} = {total_removals}")

print(f"\nThus, the minimum number of chess pieces that must be removed is {total_removals}.")

# Example of a board with 14 removals (removing row 4 and column 5)
# X = piece, O = removed
board = [['X'] * COLS for _ in range(ROWS)]
# Remove row 4 (index 3)
for c in range(COLS):
    board[3][c] = 'O'
# Remove column 5 (index 4)
for r in range(ROWS):
    board[r][4] = 'O'

print("\nExample configuration of remaining pieces (X) and removed pieces (O):")
for r in range(ROWS):
    print(" ".join(board[r]))
