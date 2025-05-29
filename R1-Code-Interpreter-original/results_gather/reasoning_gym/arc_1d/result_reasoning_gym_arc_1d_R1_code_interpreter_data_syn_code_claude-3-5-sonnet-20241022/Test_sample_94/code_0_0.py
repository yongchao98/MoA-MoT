def process_grid(input_str):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_str.split()]
    
    # Find first non-zero number position
    first_non_zero = 0
    while first_non_zero < len(numbers) and numbers[first_non_zero] == 0:
        first_non_zero += 1
    
    # Find last non-zero number position
    last_non_zero = len(numbers) - 1
    while last_non_zero >= 0 and numbers[last_non_zero] == 0:
        last_non_zero -= 1
    
    # Extract non-zero sequence
    non_zero_sequence = numbers[first_non_zero:last_non_zero + 1]
    
    # Create output with 4 more leading zeros
    output = [0] * (first_non_zero + 4)  # Original leading zeros + 4 more
    output.extend(non_zero_sequence)
    output.extend([0] * (len(numbers) - len(output)))  # Fill remaining with zeros
    
    # Convert to string
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 8 3 3 1 9 2 5 7 5 2 9 2 8 3 2 6 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))