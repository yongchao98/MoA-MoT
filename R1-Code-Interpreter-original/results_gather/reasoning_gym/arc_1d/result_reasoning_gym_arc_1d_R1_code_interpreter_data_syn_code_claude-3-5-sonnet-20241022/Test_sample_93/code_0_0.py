def analyze_transformation(input_grid, output_grid):
    # Find non-zero sequences and their positions
    def find_sequences(grid):
        sequences = []
        current_seq = []
        start_pos = None
        
        for i, num in enumerate(grid):
            if num != 0:
                if not current_seq:
                    start_pos = i
                current_seq.append(num)
            elif current_seq:
                sequences.append((current_seq.copy(), start_pos))
                current_seq = []
                start_pos = None
        
        if current_seq:
            sequences.append((current_seq.copy(), start_pos))
        return sequences

    # Analyze input and output sequences
    input_seqs = find_sequences(input_grid)
    output_seqs = find_sequences(output_grid)
    
    # Print analysis
    print("Input sequences:", input_seqs)
    print("Output sequences:", output_seqs)
    
    # Test input
    test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 7, 0]
    
    # Apply transformation rule:
    # 1. Find non-zero sequences
    # 2. Rotate them according to pattern
    test_sequences = find_sequences(test_input)
    
    # Create output based on pattern
    result = [0] * len(test_input)
    if test_sequences:
        seq, pos = test_sequences[0]
        # Place sequence earlier in the output
        new_pos = 7  # Based on observed pattern
        for i, num in enumerate(seq):
            result[new_pos + i] = num
    
    print("Predicted output:", result)

# Example 1
input1 = [0, 0, 2, 7, 9, 8, 3, 2, 1, 5, 3, 6, 0, 0, 0, 0, 0]
output1 = [9, 8, 3, 2, 1, 5, 3, 6, 0, 0, 0, 0, 0, 0, 0, 2, 7]

analyze_transformation(input1, output1)