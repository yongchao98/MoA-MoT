def transform_grid(input_grid):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_grid.split()]
    
    # Get non-zero number and its positions
    non_zeros = [x for x in numbers if x != 0]
    first_num = non_zeros[0]
    
    # Create the output pattern
    output = [first_num]  # First number
    output.append(0)      # Zero in second position
    
    # Add two non-zero numbers
    output.extend([non_zeros[1], non_zeros[2]])
    
    # Add a zero
    output.append(0)
    
    # Add three non-zero numbers
    output.extend(non_zeros[3:6])
    
    # Fill remaining positions with zeros
    while len(output) < len(numbers):
        output.append(0)
    
    return ' '.join(map(str, output))

# Test input
test_input = "9 0 9 9 9 0 9 9 0 0 0 0 0"
print(transform_grid(test_input))