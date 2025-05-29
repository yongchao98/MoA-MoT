def analyze_pattern(input_grid, output_grid):
    # Convert string to list of integers
    input_counts = [input_grid.count(i) for i in range(3)]
    output_counts = [output_grid.count(i) for i in range(3)]
    
    # Print input and output analysis
    print("Input counts (0,1,2):", input_counts)
    print("Output counts (0,1,2):", output_counts)
    print("Output is sorted:", output_grid == sorted(output_grid))

# Test examples
example1_in = [1,1,1,1,2,1,0,0,1,2,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0]
example1_out = [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2]

example2_in = [0,0,0,2,2,2,2,1,2,0,1,1,0,0,0,2,1,0,0,2,2,2,1,1,2,2]
example2_out = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2]

print("Example 1:")
analyze_pattern(example1_in, example1_out)
print("\nExample 2:")
analyze_pattern(example2_in, example2_out)

# Now analyze test input
test_input = [0,0,0,2,0,2,1,0,0,0,1,1,1,0,2,1,0,0,1,1,0,2,0,0,0,0]
counts = [test_input.count(i) for i in range(3)]
print("\nTest input counts (0,1,2):", counts)

# Generate output based on counts
output = [0] * counts[0] + [1] * counts[1] + [2] * counts[2]
print("Predicted output:", output)