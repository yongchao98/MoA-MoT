def analyze_pattern(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    def find_sequences(arr):
        sequences = []
        current_seq = []
        current_val = None
        for i, val in enumerate(arr):
            if val != 0:
                if current_val is None or val == current_val:
                    current_seq.append(i)
                    current_val = val
                else:
                    if current_seq:
                        sequences.append((current_val, current_seq))
                    current_seq = [i]
                    current_val = val
            else:
                if current_seq:
                    sequences.append((current_val, current_seq))
                    current_seq = []
                    current_val = None
        if current_seq:
            sequences.append((current_val, current_seq))
        return sequences

    input_seqs = find_sequences(input_arr)
    output_seqs = find_sequences(output_arr)
    
    print(f"Input sequences: {input_seqs}")
    print(f"Output sequences: {output_seqs}")

# Test with examples
print("Example 1:")
analyze_pattern("0 0 0 0 0 0 0 0 6 6 6 6 0 0 0 0", 
                "0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 0")

print("\nExample 2:")
analyze_pattern("0 0 0 6 6 6 6 6 6 6 6 6 6 6 0 0",
                "6 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6")

print("\nExample 3:")
analyze_pattern("0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0",
                "0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 0")

# Analyze test input
test_input = "9 9 0 0 0 0 0 0 0 0 0 0 0 0 0 9"
print("\nTest input:")
analyze_pattern(test_input, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")