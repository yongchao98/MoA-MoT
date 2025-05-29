def find_pattern(input_grid):
    # Convert string to list of integers
    numbers = [int(x) for x in input_grid.split()]
    
    # Count leading zeros
    leading_zeros = 0
    for num in numbers:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Get non-zero sequence
    non_zeros = [x for x in numbers if x != 0]
    
    # Create output with one more leading zero
    output = [0] * (leading_zeros + 1)
    output.extend(non_zeros)
    # Pad with zeros if needed
    output.extend([0] * (19 - len(output)))
    
    return output[:19]  # Ensure length is 19

# Test input
test_input = "0 9 2 1 2 8 6 6 9 8 0 0 0 0 0 0 0 0 0"
result = find_pattern(test_input)
print(" ".join(str(x) for x in result))