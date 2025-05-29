def transform_grid(input_grid):
    output_grid = []
    i = 0
    while i < len(input_grid):
        if input_grid[i] != 0:
            # Start of a non-zero sequence
            while i < len(input_grid) and input_grid[i] != 0:
                output_grid.append(input_grid[i])
                i += 1
            # Insert '3' after the sequence
            output_grid.append(3)
        else:
            # Zero element
            output_grid.append(input_grid[i])
            i += 1
    return output_grid

# Test input grid
input_grid = [0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 5, 5]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")