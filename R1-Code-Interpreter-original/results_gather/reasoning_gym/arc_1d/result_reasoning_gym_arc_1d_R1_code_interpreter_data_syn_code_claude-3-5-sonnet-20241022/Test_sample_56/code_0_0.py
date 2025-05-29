def analyze_and_predict(input_grid):
    # Convert input string to list
    grid = [int(x) for x in input_grid.split()]
    
    # Find sequences of non-zero numbers
    sequences = []
    current_seq = {'num': None, 'start': None, 'length': 0}
    
    for i, num in enumerate(grid):
        if num != 0:
            if current_seq['num'] == num:
                current_seq['length'] += 1
            else:
                if current_seq['num'] is not None:
                    sequences.append(current_seq.copy())
                current_seq = {'num': num, 'start': i, 'length': 1}
        elif current_seq['num'] is not None:
            sequences.append(current_seq.copy())
            current_seq = {'num': None, 'start': None, 'length': 0}
    
    if current_seq['num'] is not None:
        sequences.append(current_seq.copy())
    
    # Create output grid
    output = [0] * 16
    
    # For each sequence, shift it right by the appropriate amount
    for seq in sequences:
        num = seq['num']
        start = seq['start']
        length = seq['length']
        
        # Shift right by length
        new_start = start + length
        for i in range(length):
            if new_start + i < 16:  # Ensure we don't go out of bounds
                output[new_start + i] = num
    
    return ' '.join(str(x) for x in output)

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0"
result = analyze_and_predict(test_input)
print(result)