def transform_grid_optimized(input_grid):
    frequency = {}
    most_frequent = None
    max_count = 0
    output_grid = []

    for num in input_grid:
        if num != 0:
            # Update frequency count
            if num in frequency:
                frequency[num] += 1
            else:
                frequency[num] = 1

            # Check if this is the new most frequent number
            if frequency[num] > max_count:
                most_frequent = num
                max_count = frequency[num]

    # Build the output grid using the most frequent number
    for num in input_grid:
        if num == 0:
            output_grid.append(0)
        else:
            output_grid.append(most_frequent)

    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 9, 2, 2, 4, 1, 1, 1, 1, 4, 1, 1, 9, 1, 1, 1, 1, 0, 0]
output_grid = transform_grid_optimized(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")