def process_grid(input_grid):
    # Extract non-zero numbers
    non_zeros = [x for x in input_grid if x != 0]
    
    # Sort the non-zero numbers
    sorted_numbers = sorted(non_zeros)
    
    # Create output grid: sorted non-zeros followed by zeros
    output = sorted_numbers + [0] * (len(input_grid) - len(sorted_numbers))
    
    # Print the result
    print(' '.join(map(str, output)))

# Test input
test_input = [0, 8, 7, 0, 3, 4, 0, 3, 1, 0, 9, 0, 0, 0, 0, 9, 2, 7, 5, 4]
process_grid(test_input)