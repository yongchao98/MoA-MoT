def find_sequences(grid):
    sequences = []
    current_seq = []
    for i, num in enumerate(grid):
        if num != 0:
            current_seq.append((i, num))
        elif current_seq:
            if len(current_seq) > 1:
                sequences.append(current_seq)
            current_seq = []
    if current_seq and len(current_seq) > 1:
        sequences.append(current_seq)
    return sequences

def transform_grid(input_grid):
    output = [0] * len(input_grid)
    sequences = find_sequences(input_grid)
    
    # Preserve sequences
    for seq in sequences:
        for pos, val in seq:
            output[pos] = val
    
    # Find standalone numbers
    used_positions = set(pos for seq in sequences for pos, _ in seq)
    standalone = [(i, num) for i, num in enumerate(input_grid) 
                 if num != 0 and i not in used_positions]
    
    # Apply transformation for standalone numbers
    for pos, num in standalone:
        # Replicate in pattern of 3
        for i in range(3):
            if pos + i < len(output):
                output[pos + i] = num
    
    input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    result = transform_grid(input_grid)
    print(' '.join(map(str, result)))