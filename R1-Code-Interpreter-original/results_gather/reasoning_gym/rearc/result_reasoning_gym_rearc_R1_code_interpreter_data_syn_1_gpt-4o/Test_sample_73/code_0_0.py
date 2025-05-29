def analyze_transformation(input_grid, output_grid):
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    output_rows = len(output_grid)
    output_cols = len(output_grid[0])
    
    # Check which columns are retained
    retained_columns = []
    for col in range(input_cols):
        if any(input_grid[row][col] in [1, 4, 6, 8] for row in range(input_rows)):
            retained_columns.append(col)
    
    # Check the transformation rule
    transformation_rule = {}
    for row in range(output_rows):
        for col in range(output_cols):
            if output_grid[row][col] == 3:
                transformation_rule[(row, col)] = 3
            else:
                transformation_rule[(row, col)] = output_grid[row][col]
    
    return retained_columns, transformation_rule

# Example 1
input_grid_1 = [
    [0, 0, 0, 1, 6, 6, 6],
    [0, 0, 0, 1, 6, 6, 6],
    [0, 0, 0, 1, 6, 6, 6],
    [0, 0, 0, 1, 6, 0, 6],
    [0, 4, 0, 1, 6, 6, 6]
]
output_grid_1 = [
    [3, 3, 3],
    [3, 3, 3],
    [3, 3, 3],
    [3, 0, 3],
    [3, 0, 3]
]

# Example 2
input_grid_2 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 4, 1],
    [1, 1, 1, 1, 1],
    [8, 8, 8, 8, 8],
    [1, 1, 1, 1, 5],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1]
]
output_grid_2 = [
    [1, 1, 1, 1, 3],
    [1, 1, 1, 3, 1],
    [1, 1, 1, 1, 1]
]

# Example 3
input_grid_3 = [
    [9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9],
    [9, 9, 6, 9, 9],
    [4, 4, 4, 4, 4],
    [9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9],
    [9, 9, 9, 2, 9]
]
output_grid_3 = [
    [9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9],
    [9, 9, 3, 3, 9]
]

# Example 4
input_grid_4 = [
    [6, 6, 6, 6, 0, 1, 1, 1, 1],
    [6, 6, 1, 6, 0, 1, 1, 1, 1],
    [6, 6, 6, 6, 0, 1, 1, 1, 1],
    [6, 1, 6, 6, 0, 1, 1, 1, 1],
    [6, 6, 6, 6, 0, 2, 1, 1, 1],
    [6, 6, 6, 6, 0, 1, 1, 1, 2]
]
output_grid_4 = [
    [3, 3, 3, 3],
    [3, 3, 1, 3],
    [3, 3, 3, 3],
    [3, 1, 3, 3],
    [1, 3, 3, 3],
    [3, 3, 3, 1]
]

# Analyze transformations
retained_columns_1, transformation_rule_1 = analyze_transformation(input_grid_1, output_grid_1)
retained_columns_2, transformation_rule_2 = analyze_transformation(input_grid_2, output_grid_2)
retained_columns_3, transformation_rule_3 = analyze_transformation(input_grid_3, output_grid_3)
retained_columns_4, transformation_rule_4 = analyze_transformation(input_grid_4, output_grid_4)

print("Example 1 - Retained Columns:", retained_columns_1)
print("Example 1 - Transformation Rule:", transformation_rule_1)
print("Example 2 - Retained Columns:", retained_columns_2)
print("Example 2 - Transformation Rule:", transformation_rule_2)
print("Example 3 - Retained Columns:", retained_columns_3)
print("Example 3 - Transformation Rule:", transformation_rule_3)
print("Example 4 - Retained Columns:", retained_columns_4)
print("Example 4 - Transformation Rule:", transformation_rule_4)