def transform_grid(input_grid):
    # Step 1: Count the frequency of each element
    from collections import Counter
    frequency = Counter(input_grid)
    
    # Step 2: Sort elements by value in descending order
    sorted_elements = sorted(frequency.keys(), reverse=True)
    
    # Step 3: Construct the output grid
    output_grid = []
    for element in sorted_elements:
        output_grid.extend([element] * frequency[element])
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")