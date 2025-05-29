def transform_grid(input_grid):
    # Count occurrences of each number in the input grid
    frequency = {}
    for number in input_grid:
        if number in frequency:
            frequency[number] += 1
        else:
            frequency[number] = 1

    # Determine the order of numbers in the output grid
    # Based on the examples, we prioritize 2s, then 0s, then 1s
    output_grid = []
    if 2 in frequency:
        output_grid.extend([2] * frequency[2])
    if 0 in frequency:
        output_grid.extend([0] * (len(input_grid) - frequency[2] - frequency.get(1, 0)))
    if 1 in frequency:
        output_grid.extend([1] * frequency[1])

    return output_grid

# Test input
input_grid = [2, 2, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")