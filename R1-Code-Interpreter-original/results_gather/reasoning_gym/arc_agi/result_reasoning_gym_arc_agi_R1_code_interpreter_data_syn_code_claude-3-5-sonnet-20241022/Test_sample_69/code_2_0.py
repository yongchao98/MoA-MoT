def find_blocks(grid):
    blocks = []
    rows = len(grid)
    cols = len(grid[0])
    
    # For each row of blocks (skipping separator rows)
    for block_row in range(1, rows-2, 5):
        row_blocks = []
        # For each block in the row (skipping separator columns)
        for block_col in range(1, cols-2, 4):
            # Get the 2x2 middle part of each block
            block_values = [
                grid[block_row][block_col], grid[block_row][block_col+1],
                grid[block_row+1][block_col], grid[block_row+1][block_col+1]
            ]
            # Check if all values in the 2x2 block are the same
            if len(set(block_values)) == 1 and block_values[0] != 3:
                row_blocks.append(block_values[0])
            # If not all same, but not all 3s, find the non-3 value
            elif set(block_values) != {3}:
                non_3_values = [x for x in block_values if x != 3]
                if len(set(non_3_values)) == 1:
                    row_blocks.append(non_3_values[0])
        if row_blocks:
            blocks.append(row_blocks)

    return blocks

# Test input
test_input = [
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [3,2,2,3,5,3,2,2,3,5,3,5,5,3,5,3,5,5,3,5,3],
    [3,2,2,3,5,3,2,2,3,5,3,5,5,3,5,3,5,5,3,5,3],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [3,5,5,3,5,3,0,0,3,5,3,2,2,3,5,3,5,5,3,5,3],
    [3,5,5,3,5,3,0,0,3,5,3,2,2,3,5,3,5,5,3,5,3],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [3,0,0,3,5,3,5,5,3,5,3,5,5,3,5,3,1,1,3,5,3],
    [3,0,0,3,5,3,5,5,3,5,3,5,5,3,5,3,1,1,3,5,3],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3],
    [3,4,4,3,5,3,0,0,3,5,3,1,1,3,5,3,2,2,3,5,3],
    [3,4,4,3,5,3,0,0,3,5,3,1,1,3,5,3,2,2,3,5,3],
    [3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3,3,3,3,5,3]
]

blocks = find_blocks(test_input)
print("Blocks by row:")
for row in blocks:
    print(row)