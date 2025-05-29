def transform_grid(grid):
    result = list(grid)
    n = len(grid)
    
    # Find the first non-zero sequence
    start = 0
    while start < n and result[start] == 0:
        start += 1
    
    if start >= n:
        return result
        
    # Find end of sequence
    end = start
    while end < n and result[end] != 0:
        end += 1
    
    # Now we have the full sequence from start to end
    sequence = result[start:end]
    if len(sequence) < 2:
        return result
        
    # Find the main repeated number
    main_num = sequence[1]  # Take second number as reference
    
    # Check if first number is different
    if sequence[0] != main_num:
        # Move first number to end
        new_sequence = sequence[1:] + [sequence[0]]
        result[start:end] = new_sequence
    # Check if last number is different
    elif sequence[-1] != main_num:
        # Move last number to front
        new_sequence = [sequence[-1]] + sequence[:-1]
        result[start:end] = new_sequence
    
    return result

# Test input
test_input = [4,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
output = transform_grid(test_input)
print(' '.join(map(str, output)))