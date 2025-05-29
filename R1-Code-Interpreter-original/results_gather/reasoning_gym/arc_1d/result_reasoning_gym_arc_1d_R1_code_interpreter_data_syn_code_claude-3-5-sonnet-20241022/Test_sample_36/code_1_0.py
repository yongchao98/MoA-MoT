def find_and_transform(grid):
    output = grid.copy()
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
            # Transform exactly 5 numbers to 1's
            for j in range(start_idx, start_idx + 5):
                output[j] = 1
            break
            
    return output

# Test input
test_input = [0, 0, 0, 9, 9, 9, 9, 0, 0, 9, 9, 9, 9, 9, 9, 0, 0, 0, 9, 9, 9, 9, 0, 0, 0, 0]

# Get and print output
output = find_and_transform(test_input)
print(' '.join(map(str, output)))