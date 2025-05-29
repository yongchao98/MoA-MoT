def extract_bands(input_grid):
    # Initialize an empty list to store the bands
    bands = []
    
    # Iterate over each row in the input grid
    for row in input_grid:
        # Convert the row to a set to get unique elements, excluding 7
        unique_elements = list(set(row) - {7})
        
        # Sort the unique elements to maintain order
        unique_elements.sort()
        
        # If this unique set of elements is not empty and not already in bands, add it
        if unique_elements and unique_elements not in bands:
            bands.append(unique_elements)
    
    return bands

# Define the test input grid
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

# Extract the bands
output_grid = extract_bands(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))