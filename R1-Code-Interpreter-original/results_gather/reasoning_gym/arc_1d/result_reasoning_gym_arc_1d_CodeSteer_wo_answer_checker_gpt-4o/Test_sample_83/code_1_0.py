def transform_grid(input_grid):
    from collections import Counter
    
    # Count occurrences of each number in the input grid
    count = Counter(input_grid)
    
    # Remove zero from the count as we don't need to transform zeros
    if 0 in count:
        del count[0]
    
    # Sort the numbers based on their value (or any other criteria if needed)
    sorted_numbers = sorted(count.keys())
    
    # Construct the output grid
    output_grid = []
    for number in sorted_numbers:
        output_grid.extend([number] * count[number])
    
    # Fill the rest with zeros to maintain the original grid size
    output_grid.extend([0] * (len(input_grid) - len(output_grid)))
    
    return output_grid

# Test the function with the given input
test_input = [2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 1]
output_grid = transform_grid(test_input)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")