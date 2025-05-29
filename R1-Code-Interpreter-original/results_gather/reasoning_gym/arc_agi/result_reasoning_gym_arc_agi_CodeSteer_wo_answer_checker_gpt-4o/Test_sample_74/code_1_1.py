def find_prominent_numbers(input_grid):
    # Flatten the grid and count occurrences of each number
    from collections import Counter
    flat_grid = [num for row in input_grid for num in row]
    counts = Counter(flat_grid)
    
    # Assume the most common number is the background
    background_number = counts.most_common(1)[0][0]
    
    # Identify prominent numbers (those that are not the background)
    prominent_numbers = [num for num, count in counts.items() if num != background_number]
    
    return prominent_numbers

def find_bounding_box(input_grid, prominent_numbers):
    min_row, max_row = float('inf'), -float('inf')
    min_col, max_col = float('inf'), -float('inf')
    
    for r, row in enumerate(input_grid):
        for c, num in enumerate(row):
            if num in prominent_numbers:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    return min_row, max_row, min_col, max_col

def extract_output_grid(input_grid):
    prominent_numbers = find_prominent_numbers(input_grid)
    min_row, max_row, min_col, max_col = find_bounding_box(input_grid, prominent_numbers)
    
    # Extract the subgrid
    output_grid = [row[min_col:max_col+1] for row in input_grid[min_row:max_row+1]]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 9, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 9, 4, 9, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 9, 9, 9, 9, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 9, 9, 4, 9, 9, 9, 9, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 9, 4, 4, 9, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 6, 4, 4, 4, 6, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 6, 4, 4, 6, 4, 6, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 6, 4, 6, 4, 6, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4],
    [4, 4, 4, 6, 4, 6, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 4, 4],
    [4, 4, 6, 4, 6, 6, 4, 6, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 5, 5, 4, 4],
    [4, 4, 6, 4, 4, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
]

output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))