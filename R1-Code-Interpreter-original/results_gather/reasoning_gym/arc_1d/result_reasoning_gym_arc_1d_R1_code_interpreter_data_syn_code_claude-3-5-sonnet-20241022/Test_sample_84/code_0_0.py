def find_output_grid(input_grid):
    n = len(input_grid)
    output = input_grid.copy()
    
    # Function to check if position should be marked with 3
    def should_mark_three(pos):
        if input_grid[pos] != 0:  # If position is not 0, don't mark
            return False
            
        # Check all adjacent positions (including diagonals)
        for offset in [-1, 0, 1]:
            for i in range(max(0, pos-1), min(n, pos+2)):
                if i == pos:
                    continue
                # If adjacent to non-zero number, mark with 3
                if i >= 0 and i < n and input_grid[i] != 0:
                    return True
        return False
    
    # Process each position
    for i in range(n):
        if should_mark_three(i):
            output[i] = 3
            
    return output

# Test input
test_input = [0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 5, 5]
result = find_output_grid(test_input)
print(' '.join(map(str, result)))