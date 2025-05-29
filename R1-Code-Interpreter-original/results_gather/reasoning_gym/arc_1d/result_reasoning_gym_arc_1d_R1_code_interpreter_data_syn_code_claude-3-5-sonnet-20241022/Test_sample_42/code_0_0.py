def transform_grid(input_grid):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_grid.split()]
    
    # Get the first non-zero number
    first_num = next(x for x in numbers if x != 0)
    
    # Count non-zero numbers
    non_zeros = [x for x in numbers if x != 0]
    
    # Create output grid
    output = [first_num]  # First non-zero number
    output.append(0)      # Zero in second position
    
    # Add remaining non-zero numbers consecutively
    output.extend(non_zeros[1:])
    
    # Fill remaining positions with zeros
    while len(output) < len(numbers):
        output.append(0)
    
    return ' '.join(map(str, output))

# Test input
test_input = "9 0 9 9 9 0 9 9 0 0 0 0 0"
print(transform_grid(test_input))