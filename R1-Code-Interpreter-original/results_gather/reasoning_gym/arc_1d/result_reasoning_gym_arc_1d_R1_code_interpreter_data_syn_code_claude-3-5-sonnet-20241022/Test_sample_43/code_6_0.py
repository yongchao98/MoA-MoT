def transform_grid(input_str):
    numbers = [int(x) for x in input_str.split()]
    
    # Find first and last non-zero positions
    first_nonzero = next((i for i, x in enumerate(numbers) if x != 0), len(numbers))
    last_nonzero = len(numbers) - next((i for i, x in enumerate(reversed(numbers)) if x != 0), len(numbers)) - 1
    
    result = numbers.copy()
    if first_nonzero <= last_nonzero:
        # Extract non-zero sequence
        sequence = numbers[first_nonzero:last_nonzero+1]
        
        # Split into groups of 3
        groups = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        
        # Reverse each group
        reversed_groups = [group[::-1] for group in groups]
        
        # Reverse the order of groups
        reversed_groups = reversed_groups[::-1]
        
        # Flatten the groups
        transformed = [num for group in reversed_groups for num in group]
        
        # Put back into result
        result[first_nonzero:last_nonzero+1] = transformed
    
    return ' '.join(map(str, result))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 3 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))