def find_first_sequence(grid):
    current_num = 0
    count = 0
    start_idx = -1
    
    for i in range(len(grid)):
        if grid[i] == 0:
            current_num = 0
            count = 0
            continue
            
        if grid[i] == current_num:
            count += 1
        else:
            current_num = grid[i]
            count = 1
            start_idx = i
            
        if count >= 5:
            return (start_idx, count, current_num)
            
    return None

# Test input
test_input = [0, 0, 0, 9, 9, 9, 9, 0, 0, 9, 9, 9, 9, 9, 9, 0, 0, 0, 9, 9, 9, 9, 0, 0, 0, 0]

# Create output by copying input
output = test_input.copy()

# Find and transform the first sequence
sequence = find_first_sequence(test_input)
if sequence:
    start, length, num = sequence
    for i in range(start, start + length):
        output[i] = 1

# Print the output
print(' '.join(map(str, output)))