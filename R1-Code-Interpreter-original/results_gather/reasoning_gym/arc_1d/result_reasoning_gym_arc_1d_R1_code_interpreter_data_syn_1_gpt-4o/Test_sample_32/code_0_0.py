def transform_grid(input_grid):
    output_grid = []
    first_block_converted = False
    for num in input_grid:
        if num == 3 and not first_block_converted:
            output_grid.append(2)
        elif num == 3 and first_block_converted:
            output_grid.append(1)
        else:
            output_grid.append(num)
        
        if num == 0:
            first_block_converted = True
    
    return output_grid

# Test input
input_grid = [3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)