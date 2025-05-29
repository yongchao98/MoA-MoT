def analyze_pattern(input_grid, output_grid):
    # Convert string grids to 2D lists of integers
    input_rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    output_rows = [[int(x) for x in row.split()] for row in output_grid.strip().split('\n')]
    
    # Find which columns from input match the output
    matches = []
    for col in range(len(input_rows[0])-2):
        match = True
        for row in range(len(input_rows)):
            if [input_rows[row][col], input_rows[row][col+1], input_rows[row][col+2]] != output_rows[row]:
                match = False
                break
        if match:
            matches.append(col)
    return matches

# Example 1
input1 = """6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 9 6 5 6 6 6 5 6 6 6 5 6 6 6
6 9 6 5 6 9 6 5 6 6 6 5 6 6 6
6 4 6 5 6 9 6 5 6 9 6 5 6 6 6
6 6 6 5 6 4 6 5 6 9 6 5 6 9 6
6 6 6 5 6 6 6 5 6 4 6 5 6 9 6
6 6 6 5 6 6 6 5 6 6 6 5 6 4 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6
6 6 6 5 6 6 6 5 6 6 6 5 6 6 6"""

output1 = """6 6 6
6 6 6
6 6 6
6 6 6
6 6 6
6 9 6
6 9 6
6 4 6
6 6 6
6 6 6
6 6 6
6 6 6
6 6 6"""

matches = analyze_pattern(input1, output1)
print(f"Matching columns start at positions: {matches}")

# Now let's process the test input and generate output
test_input = """6 6 4 5 6 6 4 5 6 6 4 5 6 6 4
6 4 4 5 6 4 4 5 6 4 4 5 6 4 4
9 9 9 5 4 4 4 5 4 4 4 5 4 4 4
6 9 9 5 6 9 9 5 6 4 4 5 6 4 4
6 6 9 5 6 6 9 5 6 6 9 5 6 6 4
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
9 9 9 5 9 9 9 5 9 9 9 5 9 9 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
6 6 9 5 6 6 9 5 6 6 9 5 6 6 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
9 9 9 5 9 9 9 5 9 9 9 5 9 9 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
6 6 9 5 6 6 9 5 6 6 9 5 6 6 9"""

# Extract the pattern (first three columns) and print the result
input_rows = [[int(x) for x in row.split()] for row in test_input.strip().split('\n')]
result = '\n'.join(' '.join(str(x) for x in row[:3]) for row in input_rows)
print("\nPredicted output:")
print(result)