def extract_bands(input_grid):
    bands = []
    current_band = set()
    previous_row = None

    for row in input_grid:
        # Convert row to a set of unique numbers
        row_set = set(row)
        
        # If the row set is different from the previous row set, it's a new band
        if row_set != previous_row:
            if current_band:
                bands.append(list(current_band))
            current_band = row_set
        previous_row = row_set

    # Add the last band
    if current_band:
        bands.append(list(current_band))

    return bands

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 5, 5, 5, 5, 7, 7],
    [7, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Extract bands and print the output grid
output_grid = extract_bands(input_grid)
for row in output_grid:
    print(" ".join(map(str, row)))