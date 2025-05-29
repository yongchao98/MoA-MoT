def transform_grid(input_str):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_str.split()]
    
    # Separate non-zero numbers and count zeros
    non_zeros = [x for x in numbers if x != 0]
    zero_count = numbers.count(0)
    
    # Create output by combining non-zeros and zeros
    output = non_zeros + [0] * zero_count
    
    # Print result
    print(' '.join(map(str, output)))

# Test input
test_input = "0 0 0 0 0 1 5 0 0 0"
transform_grid(test_input)