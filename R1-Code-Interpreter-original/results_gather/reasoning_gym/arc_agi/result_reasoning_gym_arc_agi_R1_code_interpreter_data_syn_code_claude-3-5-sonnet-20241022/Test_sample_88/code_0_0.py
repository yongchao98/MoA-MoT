def find_vertical_sequences(grid):
    height = len(grid)
    width = len(grid[0])
    sequences = []
    
    # For each column
    for col in range(width):
        # Find sequences of 3s and 6s
        current_num = None
        sequence_length = 0
        sequence_start = 0
        
        for row in range(height):
            if grid[row][col] in [3, 6]:
                if current_num is None:
                    current_num = grid[row][col]
                    sequence_length = 1
                    sequence_start = row
                elif grid[row][col] == current_num:
                    sequence_length += 1
                    if sequence_length >= 3:  # We found a valid sequence
                        sequences.append((col, current_num))
                        break
                else:
                    current_num = grid[row][col]
                    sequence_length = 1
                    sequence_start = row
            else:
                current_num = None
                sequence_length = 0
    
    # Sort sequences by column and number
    sequences = sorted(list(set(sequences)))
    return sequences

# Parse input grid
input_grid = [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,6,0,0,0,0,0,0,0],
    [0,3,0,0,0,0,0,6,0,0,0,3,0,0,0],
    [0,3,0,0,0,3,0,6,0,0,0,3,0,0,0],
    [0,3,0,0,0,3,0,6,0,0,0,3,0,0,0],
    [0,3,0,0,0,0,0,0,0,0,0,3,0,6,0],
    [0,0,0,0,0,0,0,0,0,0,0,3,0,6,0],
    [0,0,0,6,0,0,0,0,0,0,0,3,0,6,0],
    [0,0,0,6,0,0,0,0,0,0,0,0,0,6,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]

sequences = find_vertical_sequences(input_grid)
for col, num in sequences:
    print(f"{col} {num}")