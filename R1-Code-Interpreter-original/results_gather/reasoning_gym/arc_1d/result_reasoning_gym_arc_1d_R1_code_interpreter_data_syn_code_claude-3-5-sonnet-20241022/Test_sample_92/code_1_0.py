def transform_grid(input_grid):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_grid.split()]
    
    # Process in groups of 3, removing third number
    result = []
    for i in range(0, len(numbers), 3):
        group = numbers[i:i+3]
        result.extend(group[:2])  # Take first two numbers of each group
    
    # Add trailing zero
    result.append(0)
    
    return ' '.join(map(str, result))

# Test examples to verify the rule
example1_in = "1 6 3 5 0 8 0 4 6 7 6 0 0 8 0 0 0 0 9 2 0 3 0 4 0 9 8"
example1_out = "1 6 5 0 0 4 7 6 0 8 0 0 9 2 3 0 0 9 0"
example2_in = "0 0 8 0 9 0 0 9 0 0 8 0 8 0 0 7 0 6 0 7 8 0 0 2 0 5 2"
example2_out = "0 0 0 9 0 9 0 8 8 0 7 0 0 7 0 0 0 5 0"
example3_in = "9 0 0 4 0 1 2 0 0 3 7 8 0 0 0 0 2 7 6 9 4 0 0 9 0 0 0"
example3_out = "9 0 4 0 2 0 3 7 0 0 0 2 6 9 0 0 0 0 0"

# Verify rule works for examples
print("Example 1 matches:", transform_grid(example1_in) == example1_out)
print("Example 2 matches:", transform_grid(example2_in) == example2_out)
print("Example 3 matches:", transform_grid(example3_in) == example3_out)

# Apply to test input
test_input = "8 0 0 0 0 0 0 2 0 0 0 5 2 0 0 0 3 0 8 5 0 0 0 0 0 2 8"
print("\nTest output:")
print(transform_grid(test_input))