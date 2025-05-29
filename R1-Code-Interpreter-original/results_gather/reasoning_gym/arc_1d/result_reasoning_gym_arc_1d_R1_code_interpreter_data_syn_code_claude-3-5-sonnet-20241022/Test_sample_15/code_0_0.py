def analyze_pattern(input_grid, output_grid):
    print("Analyzing patterns...")
    
    # Find sequences of same numbers
    def find_sequences(grid):
        sequences = []
        current_num = grid[0]
        current_seq = [0]
        for i in range(1, len(grid)):
            if grid[i] == current_num:
                current_seq.append(i)
            else:
                if len(current_seq) > 1:
                    sequences.append((current_num, current_seq[:]))
                current_num = grid[i]
                current_seq = [i]
        if len(current_seq) > 1:
            sequences.append((current_num, current_seq[:]))
        return sequences
    
    # Test cases
    test_cases = [
        ([0,0,0,0,2,0,7,7,7,7,7,7], [0,0,0,0,2,7,7,7,7,7,7,0]),
        ([2,0,6,6,6,6,6,6,6,6,6,6], [2,6,6,6,6,6,6,6,6,6,6,0]),
        ([8,2,0,0,0,0,0,0,0,0,0,0], [8,2,0,0,0,0,0,0,0,0,0,0])
    ]
    
    # Analyze each test case
    for input_seq, output_seq in test_cases:
        print("\nInput:", input_seq)
        print("Output:", output_seq)
        
        # Find sequences in input
        input_sequences = find_sequences(input_seq)
        print("Input sequences:", input_sequences)
        
        # Check changes
        changes = [(i, input_seq[i], output_seq[i]) for i in range(len(input_seq)) if input_seq[i] != output_seq[i]]
        print("Changes:", changes)

# Test with example cases
analyze_pattern(None, None)