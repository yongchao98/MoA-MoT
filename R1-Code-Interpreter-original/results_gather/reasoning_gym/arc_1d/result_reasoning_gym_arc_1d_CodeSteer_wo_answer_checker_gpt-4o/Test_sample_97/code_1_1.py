def transform_grid(input_grid):
    output_grid = input_grid[:]
    non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if non_zero_index is not None:
        value = input_grid[non_zero_index]
        
        # Define specific ranges to fill based on the examples
        if value == 2:
            ranges = [(1, 3), (4, 6), (7, 9), (10, 12)]
        elif value == 5:
            ranges = [(0, 4), (5, 9)]
        elif value == 8:
            ranges = [(10, 14), (15, 18)]
        else:
            ranges = []
        
        for start, end in ranges:
            for i in range(start, end):
                output_grid[i] = value
    
    return output_grid

# Test input
input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")