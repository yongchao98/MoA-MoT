def generate_output_grid(input_grid):
    # Flatten the input grid to count occurrences of each number
    from collections import Counter
    flat_input = [num for row in input_grid for num in row]
    num_counts = Counter(flat_input)
    
    # Identify the most frequent number
    most_frequent_num = num_counts.most_common(1)[0][0]
    
    # Determine the alternating number
    # Based on examples, it seems to be 2 or 5, let's choose 2 for simplicity
    alternating_num = 2
    
    # Generate the output grid
    output_grid = []
    for i in range(len(input_grid)):
        if i % 2 == 0:
            # Row with the most frequent number
            output_grid.append([most_frequent_num] * len(input_grid[0]))
        else:
            # Row with the alternating number
            output_grid.append([alternating_num] * len(input_grid[0]))
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 3, 2, 7, 7, 1, 7, 9, 6],
    [6, 1, 4, 5, 4, 6, 3, 4, 9, 6],
    [1, 4, 6, 6, 4, 8, 8, 9, 4, 3],
    [7, 3, 9, 8, 0, 0, 0, 0, 0, 4],
    [4, 4, 6, 6, 0, 0, 0, 0, 0, 2],
    [7, 3, 8, 3, 0, 0, 0, 0, 0, 5],
    [4, 4, 8, 8, 0, 0, 0, 0, 0, 9],
    [6, 3, 1, 2, 0, 1, 0, 0, 0, 3],
    [5, 7, 0, 4, 0, 0, 0, 0, 0, 6],
    [7, 4, 2, 8, 4, 7, 8, 2, 1, 9]
]

# Generate and print the output grid
output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))