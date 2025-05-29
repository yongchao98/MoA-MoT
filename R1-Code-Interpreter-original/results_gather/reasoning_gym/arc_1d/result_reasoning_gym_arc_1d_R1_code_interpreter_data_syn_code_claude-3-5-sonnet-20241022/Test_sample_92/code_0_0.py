def transform_grid(input_grid):
    # Convert input string to list of numbers
    numbers = input_grid.split()
    numbers = [int(x) for x in numbers]
    
    # Shift all numbers one position left
    output = numbers[1:] + [0]
    
    # Convert back to string format
    return ' '.join(map(str, output))

# Test examples to verify the rule
example1_in = "1 6 3 5 0 8 0 4 6 7 6 0 0 8 0 0 0 0 9 2 0 3 0 4 0 9 8"
example1_out = "1 6 3 5 8 0 4 6 7 6 0 0 8 0 0 0 0 9 2 0 3 0 4 0 9 8 0"
example2_in = "0 0 8 0 9 0 0 9 0 0 8 0 8 0 0 7 0 6 0 7 8 0 0 2 0 5 2"
example2_out = "0 8 0 9 0 0 9 0 0 8 0 8 0 0 7 0 6 0 7 8 0 0 2 0 5 2 0"

# Verify rule works for examples
print("Example 1 matches:", transform_grid(example1_in) == example1_out)
print("Example 2 matches:", transform_grid(example2_in) == example2_out)

# Apply to test input
test_input = "8 0 0 0 0 0 0 2 0 0 0 5 2 0 0 0 3 0 8 5 0 0 0 0 0 2 8"
print("\nTest output:")
print(transform_grid(test_input))