# Initial matrix
matrix = [
    [7, 0, 2],
    [4, 0, 9],
    [0, 7, 6],
    [7, 8, 4],
    [2, 9, 5]
]

# Step 2: Remove every 2-nd column (1-indexed)
matrix = [[row[0], row[2]] for row in matrix]

# Step 3: Horizontally mirror the matrix
matrix = [row[::-1] for row in matrix]

# Step 4: Set all elements divisible by 6 to zero
matrix = [[0 if element % 6 == 0 else element for element in row] for row in matrix]

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))