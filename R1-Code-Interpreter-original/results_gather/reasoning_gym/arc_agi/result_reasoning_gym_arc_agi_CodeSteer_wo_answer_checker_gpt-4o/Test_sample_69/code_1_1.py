def extract_output_grid(input_grid):
    # Flatten the grid and count occurrences of each number
    from collections import Counter
    
    flat_grid = [num for row in input_grid for num in row]
    num_counts = Counter(flat_grid)
    
    # Identify the repeating border number (most frequent)
    border_number = max(num_counts, key=num_counts.get)
    
    # Extract unique numbers excluding the border number
    unique_numbers = [num for num in flat_grid if num != border_number]
    
    # Determine the size of the output grid
    output_size = int(len(unique_numbers) ** 0.5)
    
    # Form the output grid
    output_grid = []
    for i in range(output_size):
        output_grid.append(unique_numbers[i * output_size:(i + 1) * output_size])
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3]
]

# Extract the output grid
output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))