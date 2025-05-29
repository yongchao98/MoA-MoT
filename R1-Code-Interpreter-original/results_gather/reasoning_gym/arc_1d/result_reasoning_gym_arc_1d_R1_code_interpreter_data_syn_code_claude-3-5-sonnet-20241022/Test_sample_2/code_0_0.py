def find_pattern(input_grid, output_grid):
    changes = []
    for i in range(len(input_grid)):
        if input_grid[i] != output_grid[i]:
            changes.append((i, input_grid[i], output_grid[i]))
    return changes

def solve_grid(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Create output grid
    output = grid.copy()
    
    # Find first non-zero number and its position
    for i in range(len(grid)):
        if grid[i] != 0:
            target_num = grid[i]
            # Check left side
            left_zeros = 0
            for j in range(i-1, -1, -1):
                if grid[j] == 0:
                    left_zeros += 1
                else:
                    break
            
            # Check right side
            right_zeros = 0
            next_diff_num_pos = len(grid)
            for j in range(i+1, len(grid)):
                if grid[j] == 0:
                    right_zeros += 1
                elif grid[j] != target_num:
                    next_diff_num_pos = j
                    break
            
            # Fill zeros in the direction with more consecutive zeros
            if left_zeros > right_zeros:
                for j in range(i-1, i-left_zeros-1, -1):
                    output[j] = target_num
            else:
                for j in range(i+1, next_diff_num_pos):
                    if grid[j] == 0:
                        output[j] = target_num
            
    # Print result
    print(' '.join(map(str, output)))

# Test with the given input
test_input = "0 0 3 3 3 3 3 3 3 3 3 3 3 3 2 0 0 0 0 0 0"
solve_grid(test_input)