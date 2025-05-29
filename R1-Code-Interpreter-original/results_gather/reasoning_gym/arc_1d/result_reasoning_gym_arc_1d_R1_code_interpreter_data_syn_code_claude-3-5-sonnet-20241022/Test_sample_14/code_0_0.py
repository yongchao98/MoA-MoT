def analyze_and_transform(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find non-zero sequences
    non_zero_sequences = []
    current_seq = []
    for i, num in enumerate(grid):
        if num != 0:
            current_seq.append((num, i))
        elif current_seq:
            non_zero_sequences.append(current_seq)
            current_seq = []
    if current_seq:
        non_zero_sequences.append(current_seq)
    
    # Create output grid
    output = [0] * len(grid)
    
    # Move each sequence 2 positions left
    for sequence in non_zero_sequences:
        for num, pos in sequence:
            new_pos = max(0, pos - 2)  # Ensure we don't go negative
            output[new_pos] = num
            
    # Print result
    print(" ".join(map(str, output)))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1"
analyze_and_transform(test_input)