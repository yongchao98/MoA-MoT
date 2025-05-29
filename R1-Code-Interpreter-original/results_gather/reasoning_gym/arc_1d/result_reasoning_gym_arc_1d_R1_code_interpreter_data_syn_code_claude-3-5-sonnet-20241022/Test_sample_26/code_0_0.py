def find_sequences(grid):
    sequences = []
    current_seq = [0, grid[0]]  # [start_pos, value]
    
    for i in range(1, len(grid)):
        if grid[i] != current_seq[1]:
            sequences.append([current_seq[0], i-1, current_seq[1]])  # [start, end, value]
            current_seq = [i, grid[i]]
    sequences.append([current_seq[0], len(grid)-1, current_seq[1]])
    return sequences

def shift_sequence(grid):
    sequences = find_sequences(grid)
    result = grid.copy()
    
    # Find the sequence to shift (either all same numbers or all zeros)
    for seq in sequences:
        start, end, value = seq
        length = end - start + 1
        if length > 1:  # Only consider sequences longer than 1
            # If sequence is all zeros or all same numbers
            if value == 0 or all(grid[i] == value for i in range(start, end+1)):
                # Shift left by 2
                new_start = max(0, start - 2)
                for i in range(start, end+1):
                    result[i] = grid[start]  # Reset to original value
                for i in range(new_start, new_start + length):
                    result[i] = value
                break
    
    print(' '.join(map(str, result)))

# Test input
test_input = [2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 2, 2, 2]
shift_sequence(test_input)