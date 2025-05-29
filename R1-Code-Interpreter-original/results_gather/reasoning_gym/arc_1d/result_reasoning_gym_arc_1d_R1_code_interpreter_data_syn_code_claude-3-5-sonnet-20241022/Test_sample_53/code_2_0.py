def solve_grid(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    
    # Find the majority number
    from collections import Counter
    majority = Counter(grid).most_common(1)[0][0]
    
    # Find all sequences of non-majority numbers
    sequences = []
    current_seq = []
    for i, num in enumerate(grid):
        if num != majority:
            current_seq.append((i, num))
        elif current_seq:
            sequences.append(current_seq)
            current_seq = []
    if current_seq:
        sequences.append(current_seq)
    
    # Create output grid starting with majority numbers
    output = [majority] * len(grid)
    
    # Process each sequence
    for sequence in sequences:
        # Get original positions and values
        orig_start = sequence[0][0]
        # Move sequence 2 positions left
        new_start = max(0, orig_start - 2)
        for i, (_, num) in enumerate(sequence):
            if new_start + i < len(grid):
                output[new_start + i] = num
    
    return ' '.join(map(str, output))

# Test with the given input
test_input = "0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0"
print(solve_grid(test_input))