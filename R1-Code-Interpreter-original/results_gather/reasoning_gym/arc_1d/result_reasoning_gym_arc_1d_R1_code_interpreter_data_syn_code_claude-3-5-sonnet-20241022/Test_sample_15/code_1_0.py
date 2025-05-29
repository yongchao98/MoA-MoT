def apply_rule(input_grid):
    # Find the last sequence
    last_num = input_grid[-1]
    seq_start = len(input_grid) - 1
    
    # Find where the last sequence starts
    while seq_start > 0 and input_grid[seq_start - 1] == last_num:
        seq_start -= 1
    
    # If we have a sequence of non-zero numbers at the end
    if last_num != 0 and seq_start < len(input_grid):
        # Check if there's a zero before the sequence
        if seq_start > 0 and input_grid[seq_start - 1] == 0:
            result = input_grid.copy()
            result[seq_start - 1] = last_num  # Change 0 to sequence number
            result[-1] = 0  # Change last number to 0
            return result
    
    return input_grid.copy()

# Test with the new input
test_input = [9, 9, 9, 9, 9, 9, 9, 9, 0, 0, 0, 2]
result = apply_rule(test_input)
print("Input:", test_input)
print("Output:", result)