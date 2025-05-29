def analyze_patterns():
    # Example inputs and outputs
    ex1_in = [8,8,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
    ex1_out = [0,0,0,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,0,0]
    ex2_in = [6,6,6,6,6,6,6,6,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
    ex2_out = [6,6,6,6,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
    ex3_in = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    ex3_out = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    
    test_in = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,6,0,0,0,0,0]
    
    # Analyze differences
    def find_changes(inp, out):
        return [(i, v1, v2) for i, (v1, v2) in enumerate(zip(inp, out)) if v1 != v2]
    
    print("Changes in Example 1:", find_changes(ex1_in, ex1_out))
    print("Changes in Example 2:", find_changes(ex2_in, ex2_out))
    print("Changes in Example 3:", find_changes(ex3_in, ex3_out))
    
    # Find sequences of same numbers
    def find_sequences(arr):
        sequences = []
        current_num = arr[0]
        current_count = 1
        start_pos = 0
        
        for i in range(1, len(arr)):
            if arr[i] == current_num:
                current_count += 1
            else:
                sequences.append((start_pos, current_num, current_count))
                current_num = arr[i]
                current_count = 1
                start_pos = i
        sequences.append((start_pos, current_num, current_count))
        return sequences
    
    print("\nSequences in test input:", find_sequences(test_in))

analyze_patterns()