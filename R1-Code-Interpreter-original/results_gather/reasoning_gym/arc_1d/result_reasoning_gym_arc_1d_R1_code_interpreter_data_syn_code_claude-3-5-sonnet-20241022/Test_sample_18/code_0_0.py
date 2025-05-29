def find_pattern(input_str):
    # Convert input string to list of numbers
    numbers = [int(x) for x in input_str.split()]
    
    # Create output list, initially same as input
    output = numbers.copy()
    
    # Iterate through the list
    i = 0
    while i < len(numbers):
        # If we find a 2, check if it starts a sequence of at least 2 twos
        if numbers[i] == 2:
            # Count consecutive 2s
            count = 0
            j = i
            while j < len(numbers) and numbers[j] == 2:
                count += 1
                j += 1
            
            # If we found at least 2 consecutive 2s
            if count >= 2:
                # Change first two 2s to 8s
                output[i] = 8
                output[i+1] = 8
            
            # Skip the sequence we just processed
            i = j
        else:
            i += 1
    
    # Convert back to string
    return ' '.join(map(str, output))

# Test input
test_input = "2 2 2 2 2 0 0 2 2 2 0 0 0 0 0 0 0 0"
print(find_pattern(test_input))