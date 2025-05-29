def shift_right(input_grid):
    # Convert string to list of integers
    numbers = [int(x) for x in input_grid.split()]
    
    # Get non-zero numbers
    non_zeros = [x for x in numbers if x != 0]
    
    # Create output grid
    output = [0] * len(numbers)
    
    # Place non-zero numbers shifted one position right
    for i, num in enumerate(non_zeros):
        if i + 1 < len(output):
            output[i + 1] = num
    
    return output

# Test cases
test_inputs = [
    "0 0 0 0 0 0 0 6 2 8 8 1 0 0 0 0 0 0 0",
    "0 6 9 7 7 3 1 2 2 7 3 2 3 9 8 3 7 9 0",
    "0 0 0 0 0 0 0 0 0 3 7 2 1 1 3 1 3 5 0"
]

test_outputs = [
    "0 0 0 0 0 0 0 0 6 2 8 8 1 0 0 0 0 0 0",
    "0 0 6 9 7 7 3 1 2 2 7 3 2 3 9 8 3 7 9",
    "0 0 0 0 0 0 0 0 0 0 3 7 2 1 1 3 1 3 5"
]

# Verify pattern
for i in range(len(test_inputs)):
    result = shift_right(test_inputs[i])
    expected = [int(x) for x in test_outputs[i].split()]
    print(f"Test {i+1} matches: {result == expected}")

# Apply to test input
test_input = "0 9 2 1 2 8 6 6 9 8 0 0 0 0 0 0 0 0 0"
result = shift_right(test_input)
print("\nTest input result:")
print(" ".join(str(x) for x in result))