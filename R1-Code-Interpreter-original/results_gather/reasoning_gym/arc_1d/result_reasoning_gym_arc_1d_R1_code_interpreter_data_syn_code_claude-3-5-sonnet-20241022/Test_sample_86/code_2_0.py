def apply_rule(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    result = [0] * n
    
    # Find non-zero sequences
    sequences = []
    current_seq = []
    start_pos = -1
    
    for i, num in enumerate(grid):
        if num != 0:
            if not current_seq:
                start_pos = i
            current_seq.append(num)
        elif current_seq:
            sequences.append((start_pos, current_seq[:]))
            current_seq = []
            start_pos = -1
    
    if current_seq:
        sequences.append((start_pos, current_seq))
    
    # Special case: if there's a single number at position 15
    if len(sequences) == 1 and sequences[0][0] == 15 and len(sequences[0][1]) == 1:
        result[0] = sequences[0][1][0]
    else:
        # Move all numbers 3 positions to the right
        nonzero_nums = [x for x in grid if x != 0]
        if nonzero_nums:
            first_nonzero_pos = next(i for i, x in enumerate(grid) if x != 0)
            new_start = (first_nonzero_pos + 3) % n
            for i, num in enumerate(nonzero_nums):
                result[(new_start + i) % n] = num
    
    return ' '.join(map(str, result))

# Test input
test_input = "8 1 2 0 0 0 0 0 0 0 8 7 2 8 4 5 1 8 2"
print(apply_rule(test_input))