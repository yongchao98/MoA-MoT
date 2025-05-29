def transform_grid(input_grid):
    # Create output by shifting left and adding 0 at start
    output = [0] + input_grid[:-1]
    return output

# Test all examples
example1_in = [2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0, 0]
example1_out = [0, 2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0]

example2_in = [0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 0, 3]
example2_out = [0, 0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 3]

example3_in = [5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2, 0]
example3_out = [0, 5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2]

# Verify examples
print("Example 1:", transform_grid(example1_in) == example1_out)
print("Example 2:", transform_grid(example2_in) == example2_out)
print("Example 3:", transform_grid(example3_in) == example3_out)

# Transform test input
test_input = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
result = transform_grid(test_input)
print("\nTest input:", test_input)
print("Result:", result)