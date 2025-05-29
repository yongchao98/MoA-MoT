def verify_flip_pattern(input_grid, output_grid):
    # Convert input string to 2D list
    input_rows = [[int(num) for num in row.split()] for row in input_grid.split('\n') if row]
    output_rows = [[int(num) for num in row.split()] for row in output_grid.split('\n') if row]
    
    # Check if output is vertical flip of input
    flipped_input = input_rows[::-1]
    return flipped_input == output_rows

# Test examples
example1_input = """0 3 3 0 3 3
3 0 0 3 4 0
4 4 3 0 0 6
0 4 6 4 3 3
0 3 3 6 4 4
0 6 0 3 0 4"""

example1_output = """0 6 0 3 0 4
0 3 3 6 4 4
0 4 6 4 3 3
4 4 3 0 0 6
3 0 0 3 4 0
0 3 3 0 3 3"""

example2_input = """6 6 4 0
6 0 3 4
6 6 4 3
3 3 3 3"""

example2_output = """3 3 3 3
6 6 4 3
6 0 3 4
6 6 4 0"""

# Verify pattern
print("Example 1 matches pattern:", verify_flip_pattern(example1_input, example1_output))
print("Example 2 matches pattern:", verify_flip_pattern(example2_input, example2_output))

# Generate output for test input
test_input = """4 6 6
6 4 3
0 6 3"""

# Convert test input to 2D list and flip
test_input_rows = [[int(num) for num in row.split()] for row in test_input.split('\n') if row]
test_output = test_input_rows[::-1]

# Format output
print("\nPredicted output:")
for row in test_output:
    print(' '.join(map(str, row)))