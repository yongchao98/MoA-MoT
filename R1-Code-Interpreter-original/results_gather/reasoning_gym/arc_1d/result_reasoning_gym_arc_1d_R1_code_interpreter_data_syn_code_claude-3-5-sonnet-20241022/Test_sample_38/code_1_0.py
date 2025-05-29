def find_pattern(input_str):
    # Convert input string to list
    input_list = [int(x) for x in input_str.split()]
    
    # Find first non-zero number
    first_non_zero = 0
    while first_non_zero < len(input_list) and input_list[first_non_zero] == 0:
        first_non_zero += 1
        
    # If all zeros, return all zeros
    if first_non_zero == len(input_list):
        return input_list
    
    # Find last non-zero number
    last_non_zero = len(input_list) - 1
    while last_non_zero >= 0 and input_list[last_non_zero] == 0:
        last_non_zero -= 1
    
    # Create output list
    output = [0] * len(input_list)
    
    # Copy the sequence from first non-zero to last non-zero
    for i in range(first_non_zero, last_non_zero + 1):
        output[i] = input_list[i]
    
    print(' '.join(map(str, output)))

# Test input
test_input = "0 0 0 0 6 3 5 1 7 9 1 2 4 9 6 3 2 0 0"
find_pattern(test_input)