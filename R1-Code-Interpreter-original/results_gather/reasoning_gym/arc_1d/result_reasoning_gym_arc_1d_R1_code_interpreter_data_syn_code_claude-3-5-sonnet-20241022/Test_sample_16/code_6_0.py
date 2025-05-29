def transform_grid(input_grid):
    # Find position of first non-zero number
    first_nonzero = 0
    for i, num in enumerate(input_grid):
        if num != 0:
            first_nonzero = i
            break
    
    # Create output grid
    output = []
    
    # Add zeros at start (number of zeros equals position of first non-zero + 1)
    output.extend([0] * (first_nonzero + 1))
    
    # Add remaining numbers from first non-zero to second-to-last
    non_zero_sequence = input_grid[first_nonzero:-1]
    output.extend(non_zero_sequence)
    
    # Pad with zeros if needed to match length
    while len(output) < len(input_grid):
        output.append(0)
        
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