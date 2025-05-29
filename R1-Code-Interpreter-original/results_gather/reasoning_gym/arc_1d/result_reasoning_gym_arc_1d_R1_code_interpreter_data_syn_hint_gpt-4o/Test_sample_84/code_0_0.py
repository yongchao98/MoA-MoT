def transform_grid(input_grid):
    output_grid = input_grid[:]
    last_positions = {}
    
    # Find the last occurrence of each distinct non-zero number
    for i, num in enumerate(input_grid):
        if num != 0:
            last_positions[num] = i
    
    # Insert '3' after the last occurrence of each distinct non-zero number
    offset = 0
    for num, pos in sorted(last_positions.items()):
        output_grid.insert(pos + 1 + offset, 3)
        offset += 1
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 5, 5]
output_grid = transform_grid(input_grid)
print(output_grid)