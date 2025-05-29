import numpy as np

def transform_grid(input_grid, transformation_rules):
    # Convert the input grid to a numpy array for efficient processing
    input_array = np.array(input_grid)
    
    # Create an output array initialized with the input values
    output_array = np.copy(input_array)
    
    # Apply transformation rules
    for original_value, new_value in transformation_rules.items():
        output_array[input_array == original_value] = new_value
    
    return output_array

# Define transformation rules
transformation_rules = {
    0: 0,
    6: 5,
    8: 5,
    2: 5
}

# Test input grid
input_grid = [
    [1, 7, 1, 1, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 7, 1, 1, 1],
    [1, 1, 1, 1, 1, 7, 1, 1, 7, 1, 1, 1, 7, 0, 0, 7, 0, 7],
    [1, 7, 1, 1, 7, 1, 1, 7, 7, 7, 1, 1, 0, 7, 7, 0, 0, 1],
    [1, 7, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [7, 7, 7, 7, 1, 1, 7, 1, 7, 7, 1, 7, 1, 1, 1, 7, 7, 1],
    [1, 1, 1, 1, 7, 1, 1, 7, 0, 0, 7, 0, 7, 1, 7, 7, 7, 7],
    [1, 1, 1, 1, 7, 7, 1, 7, 0, 0, 0, 0, 7, 1, 1, 1, 1, 1],
    [7, 1, 7, 7, 7, 1, 7, 7, 1, 7, 1, 1, 7, 1, 7, 7, 1, 7],
    [1, 1, 1, 1, 7, 1, 7, 7, 7, 7, 7, 1, 1, 1, 7, 1, 1, 1],
    [7, 1, 1, 1, 1, 7, 7, 1, 1, 7, 7, 1, 7, 1, 7, 1, 7, 7],
    [7, 7, 7, 7, 1, 1, 7, 1, 1, 7, 0, 7, 1, 7, 1, 7, 7, 1],
    [7, 1, 7, 7, 7, 1, 1, 1, 7, 7, 7, 0, 7, 1, 1, 7, 7, 1],
    [1, 7, 7, 7, 1, 1, 1, 7, 1, 1, 7, 7, 7, 1, 7, 7, 1, 7],
    [1, 7, 0, 0, 7, 0, 0, 7, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7],
    [7, 0, 7, 7, 0, 7, 0, 1, 1, 7, 1, 7, 7, 7, 7, 7, 7, 1],
    [1, 7, 0, 0, 7, 0, 0, 1, 7, 1, 7, 7, 1, 1, 7, 1, 7, 1],
    [1, 7, 7, 7, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7]
]

# Get the output grid
output_grid = transform_grid(input_grid, transformation_rules)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))