def transform_grid(input_grid):
    # The first element of the output is the same as the first element of the input
    first_element = input_grid[0]
    
    # The last element of the output is always 0
    last_element = 0
    
    # Find the first non-zero element after the first element
    fill_value = next((x for x in input_grid[1:] if x != 0), 0)
    
    # Construct the output grid
    output_grid = [first_element] + [fill_value] * (len(input_grid) - 2) + [last_element]
    
    return output_grid

# Test input
input_grid = [9, 9, 9, 9, 9, 9, 9, 9, 0, 0, 0, 2]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")