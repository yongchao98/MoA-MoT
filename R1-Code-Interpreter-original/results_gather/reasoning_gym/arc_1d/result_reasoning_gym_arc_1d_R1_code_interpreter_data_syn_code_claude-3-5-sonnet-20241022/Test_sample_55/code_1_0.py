def apply_rule(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    output = [0] * 16
    
    # Find sequences of non-zero numbers
    sequences = []
    current_seq = []
    current_val = None
    
    for i, val in enumerate(grid):
        if val != 0:
            if current_val is None or val == current_val:
                current_seq.append(i)
                current_val = val
            else:
                if current_seq:
                    sequences.append((current_val, current_seq[:]))
                current_seq = [i]
                current_val = val
        else:
            if current_seq:
                sequences.append((current_val, current_seq[:]))
                current_seq = []
                current_val = None
    if current_seq:
        sequences.append((current_val, current_seq))
    
    # Apply transformation rules
    for val, positions in sequences:
        # If sequence starts at beginning or end, keep it there
        if positions[0] == 0 or positions[-1] == 15:
            for pos in positions:
                output[pos] = val
        else:
            # Move sequence right by 3 positions
            for pos in positions:
                new_pos = pos + 3
                if new_pos < 16:
                    output[new_pos] = val
                else:
                    # If position would go beyond boundary, wrap to front
                    output[0] = val
    
    return ' '.join(str(x) for x in output)

# Test with all examples
print("Example 1:")
print(apply_rule("0 0 0 0 0 0 0 0 6 6 6 6 0 0 0 0"))
print("Example 2:")
print(apply_rule("0 0 0 6 6 6 6 6 6 6 6 6 6 6 0 0"))
print("Example 3:")
print(apply_rule("0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0"))
print("Test input:")
print(apply_rule("9 9 0 0 0 0 0 0 0 0 0 0 0 0 0 9"))