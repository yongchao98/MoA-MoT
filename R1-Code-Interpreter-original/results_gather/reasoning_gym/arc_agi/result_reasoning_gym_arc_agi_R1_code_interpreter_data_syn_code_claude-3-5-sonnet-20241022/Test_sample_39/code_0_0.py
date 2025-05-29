def analyze_transformations(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find where 0s change to 5s
    changes = []
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == '0' and output_grid[i][j] == '5':
                changes.append((i, j))
    
    # Analyze the pattern of changes
    print("Changes from 0 to 5 at positions (row, col):")
    for change in changes:
        print(f"Row {change[0]}, Col {change[1]}")
        
    # Print the distance between changes
    print("\nDistances between consecutive changes:")
    for i in range(1, len(changes)):
        dist = changes[i][0] - changes[i-1][0]
        print(f"Distance: {dist} rows")

# Example 1
input1 = [
    "9 0 9", "0 9 9", "9 9 0", "9 0 9", "9 9 0", "9 0 9", "9 9 0", 
    "0 9 9", "9 9 0", "9 0 9", "0 9 9", "9 9 0", "0 9 9", "9 0 9"
]
output1 = [
    "9 5 9", "0 9 9", "9 9 5", "9 0 9", "9 9 5", "9 0 9", "9 9 5",
    "0 9 9", "9 9 5", "9 0 9", "5 9 9", "9 9 0", "5 9 9", "9 0 9"
]

# Convert to 2D arrays
input_grid = [row.split() for row in input1]
output_grid = [row.split() for row in output1]

analyze_transformations(input_grid, output_grid)