def analyze_and_predict(examples, test_input):
    def find_sequence(grid):
        # Find non-zero sequences
        sequences = []
        current_seq = []
        current_num = 0
        start_pos = 0
        
        for i, num in enumerate(grid):
            if num != 0 and num == current_num:
                current_seq.append(num)
            elif num != 0:
                if current_seq:
                    sequences.append((start_pos, current_seq))
                current_seq = [num]
                current_num = num
                start_pos = i
            elif current_seq:
                sequences.append((start_pos, current_seq))
                current_seq = []
                current_num = 0
        if current_seq:
            sequences.append((start_pos, current_seq))
        return sequences

    # Analyze examples
    shifts = []
    for inp, out in examples:
        in_seq = find_sequence(inp)
        out_seq = find_sequence(out)
        if in_seq and out_seq:
            shift = in_seq[0][0] - out_seq[0][0]
            shifts.append(shift)
    
    # Determine rule
    # It appears the shift is proportional to sequence length
    test_seq = find_sequence(test_input)
    if not test_seq:
        return test_input
    
    # Create output grid
    result = [0] * 16
    seq_start, seq_nums = test_seq[0]
    seq_len = len(seq_nums)
    
    # Calculate shift based on sequence length
    shift = 1 if seq_len <= 3 else 2
    
    # Apply shift
    new_start = seq_start - shift
    for i, num in enumerate(seq_nums):
        if 0 <= new_start + i < 16:
            result[new_start + i] = num
    
    return result

# Example inputs and outputs
examples = [
    ([0,0,0,0,0,0,0,0,0,4,4,4,0,0,0,0], [0,0,0,0,0,0,0,4,4,4,0,0,0,0,0,0]),
    ([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1], [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0]),
    ([0,0,0,0,0,0,0,0,0,9,9,0,0,0,0,0], [0,0,0,0,0,0,0,9,9,0,0,0,0,0,0,0])
]

test_input = [0,0,0,0,0,0,0,0,0,0,9,9,9,9,9,9]
result = analyze_and_predict(examples, test_input)
print(' '.join(map(str, result)))