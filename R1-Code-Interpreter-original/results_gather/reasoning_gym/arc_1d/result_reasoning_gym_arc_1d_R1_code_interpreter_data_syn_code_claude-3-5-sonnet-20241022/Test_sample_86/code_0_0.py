def analyze_movement(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences and their positions
    def find_nonzero_sequences(arr):
        sequences = []
        current_seq = []
        start_pos = -1
        
        for i, num in enumerate(arr):
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
        
        return sequences
    
    input_seqs = find_nonzero_sequences(input_arr)
    output_seqs = find_nonzero_sequences(output_arr)
    
    print(f"Input sequences: {input_seqs}")
    print(f"Output sequences: {output_seqs}")

# Analyze examples
print("Example 1:")
analyze_movement("0 0 0 0 0 0 0 0 0 0 0 5 1 5 6 0 0 0 0", 
                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 1 5 6")

print("\nExample 2:")
analyze_movement("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 0 0 0",
                "9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")

print("\nExample 3:")
analyze_movement("3 9 1 4 3 7 5 2 1 9 0 0 0 0 0 0 0 0 9",
                "0 0 0 9 3 9 1 4 3 7 5 2 1 9 0 0 0 0 0")