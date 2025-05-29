def find_pattern(input_grid):
    # Convert to list to make it mutable
    output = list(input_grid)
    
    # Find pattern: 1 -> 0 -> sequence of same numbers
    for i in range(len(input_grid)-2):
        if input_grid[i] == 1 and input_grid[i+1] == 0:
            # Look for sequence of same numbers
            j = i + 2
            num = input_grid[j]
            sequence_length = 0
            
            while j < len(input_grid) and input_grid[j] == num:
                sequence_length += 1
                j += 1
            
            if sequence_length > 0:
                # Apply the transformation
                output[i+1] = num  # Replace 0 with the sequence number
                output[j-1] = 0    # Replace last number with 0
    
    return output

# Test input
test_input = [8, 8, 8, 8, 0, 0, 1, 0, 0, 0, 0]

# Test all examples to verify rule
examples = [
    [0, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 4, 4, 4, 4, 0],
    [0, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6]
]

print("Test input result:", find_pattern(test_input))
print("\nVerifying examples:")
for ex in examples:
    print(f"Input:  {ex}")
    print(f"Output: {find_pattern(ex)}")