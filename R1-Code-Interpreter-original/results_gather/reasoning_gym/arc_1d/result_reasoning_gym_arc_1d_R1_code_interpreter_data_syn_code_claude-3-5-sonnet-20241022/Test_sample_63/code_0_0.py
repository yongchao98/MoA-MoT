def transform_grid(input_grid):
    output = list(input_grid)  # Start with a copy of input
    n = len(input_grid)
    
    # For each position in the grid
    for i in range(n):
        if input_grid[i] != 0:  # Found a non-zero number
            current_value = input_grid[i]
            # Fill rightward until next non-zero or end
            j = i
            while j < n:
                if j == i or input_grid[j] == 0:  # Fill if it's the start position or a zero
                    output[j] = current_value
                else:  # Stop at next non-zero number
                    break
                j += 1
    
    return output

# Test input
test_input = [0, 0, 0, 0, 0, 9, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 7, 0, 0, 0, 5]

# Verify examples
example1_in = [6, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 5]
example1_out = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 5]

# Print result for test input
result = transform_grid(test_input)
print(" ".join(map(str, result)))

# Verify example 1 (to confirm rule is correct)
example1_result = transform_grid(example1_in)
print("\nRule verification (should match example 1):")
print("Expected:", " ".join(map(str, example1_out)))
print("Got:     ", " ".join(map(str, example1_result)))