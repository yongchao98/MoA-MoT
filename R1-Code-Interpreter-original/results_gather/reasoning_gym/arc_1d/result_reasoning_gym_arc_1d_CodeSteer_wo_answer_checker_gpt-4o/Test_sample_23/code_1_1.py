def transform_grid(input_grid):
    # Step 1: Identify the non-zero sequence
    start_index = None
    end_index = None
    for i, num in enumerate(input_grid):
        if num != 0:
            if start_index is None:
                start_index = i
            end_index = i

    # If no non-zero sequence is found, return the input as is
    if start_index is None:
        return input_grid

    # Step 2: Calculate the shift
    shift = start_index

    # Step 3: Construct the output grid
    output_grid = [0] * len(input_grid)
    sequence_length = end_index - start_index + 1
    output_grid[shift:shift + sequence_length] = input_grid[start_index:end_index + 1]

    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")