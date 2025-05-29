def transform_grid(input_grid):
    # Convert input grid to a list of integers
    input_list = list(map(int, input_grid.split()))
    
    # First column is the first element of the input grid
    first_col = input_list[0]
    
    # Find the largest and second largest unique numbers in the input grid
    unique_numbers = sorted(set(input_list), reverse=True)
    largest = unique_numbers[0]
    second_largest = unique_numbers[1] if len(unique_numbers) > 1 else 0
    
    # Construct the output grid
    output_grid = []
    output_grid.append(f"{first_col} {largest} {second_largest}")
    output_grid.append(f"{first_col} {largest} 0")
    
    # Fill the rest with zeros
    for _ in range(len(input_list) - 2):
        output_grid.append(f"{first_col} 0 0")
    
    return "\n".join(output_grid)

# Test input
test_input = """
4
4
4
4
4
4
4
4
4
4
4
4
1
4
4
4
1
9
"""

# Transform the grid and print the output
print(transform_grid(test_input))