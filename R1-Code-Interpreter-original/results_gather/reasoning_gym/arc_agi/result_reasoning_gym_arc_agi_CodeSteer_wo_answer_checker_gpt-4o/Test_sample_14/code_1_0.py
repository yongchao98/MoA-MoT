def find_pattern_center(grid):
    # Dictionary to count occurrences of each number
    number_count = {}
    
    # Scan the grid to count occurrences of each number
    for row in grid:
        for num in row:
            if num != 7:  # Assuming 7 is the background number
                if num in number_count:
                    number_count[num] += 1
                else:
                    number_count[num] = 1
    
    # Find the most frequent non-background number
    pattern_number = max(number_count, key=number_count.get)
    
    # Find the center of the pattern_number
    for i in range(1, len(grid) - 1):
        for j in range(1, len(grid[i]) - 1):
            if grid[i][j] == pattern_number:
                # Check if it's part of a pattern
                if (grid[i-1][j] == pattern_number or grid[i+1][j] == pattern_number or
                    grid[i][j-1] == pattern_number or grid[i][j+1] == pattern_number):
                    return i, j
    return None

def extract_3x3_section(grid, center):
    i, j = center
    return [row[j-1:j+2] for row in grid[i-1:i+2]]

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 5, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9],
    [7, 7, 5, 5, 5, 7, 7, 7, 7, 7, 7, 9, 9, 9],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 9, 7],
    [7, 7, 7, 7, 7, 7, 7, 8, 7, 8, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 7],
    [7, 7, 7, 9, 9, 7, 7, 7, 7, 7, 7, 5, 7, 5],
    [7, 7, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 5, 7],
    [7, 7, 7, 9, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Find the center of the pattern
center = find_pattern_center(input_grid)

# Extract the 3x3 section
output_grid = extract_3x3_section(input_grid, center)

# Print the output grid
for row in output_grid:
    print(" ".join(map(str, row)))