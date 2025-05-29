def transform_grid(input_grid):
    # Find the start and end indices of the non-zero sequence
    start_index = None
    end_index = None
    for i, value in enumerate(input_grid):
        if value != 0:
            if start_index is None:
                start_index = i
            end_index = i

    # Initialize the output grid with zeros
    output_grid = [0] * len(input_grid)

    # Calculate new start and end indices after shifting two positions to the left
    if start_index is not None and end_index is not None:
        new_start_index = max(0, start_index - 2)
        new_end_index = new_start_index + (end_index - start_index)

        # Copy the non-zero sequence to the new positions in the output grid
        for i in range(new_start_index, new_end_index + 1):
            output_grid[i] = input_grid[start_index + (i - new_start_index)]

    return output_grid

# Test input
input_grid = [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")