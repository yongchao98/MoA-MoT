def verify_pattern(input_grid, output_grid):
    # Create expected output by shifting right and adding 0 at start
    expected = [0] + input_grid[:-1]
    return expected == output_grid

# Test examples
example1_in = [2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0, 0]
example1_out = [0, 2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0]

example2_in = [0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 0, 3]
example2_out = [0, 0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 3]

example3_in = [5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2, 0]
example3_out = [0, 5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2]

# Verify all examples
print("Example 1 matches pattern:", verify_pattern(example1_in, example1_out))
print("Example 2 matches pattern:", verify_pattern(example2_in, example2_out))
print("Example 3 matches pattern:", verify_pattern(example3_in, example3_out))

# Generate output for test input
test_input = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
test_output = [0] + test_input[:-1]
print("\nPredicted output:", test_output)