def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find sequences of non-zero numbers
    def find_sequences(lst):
        sequences = []
        current_seq = []
        current_num = None
        
        for i, num in enumerate(lst):
            if num != 0:
                if num != current_num:
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = [num]
                    current_num = num
                else:
                    current_seq.append(num)
            else:
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = []
                    current_num = None
        if current_seq:
            sequences.append(current_seq)
        return sequences

    # Analyze examples
    input_seqs = find_sequences(input_list)
    output_seqs = find_sequences(output_list)
    
    print(f"Input sequences: {input_seqs}")
    print(f"Output sequences: {output_seqs}")
    print("Length comparison:", len(input_list), len(output_list))
    
# Test with example 1
print("Example 1:")
analyze_pattern("0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")

print("\nExample 2:")
analyze_pattern("4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4",
                "4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4")

print("\nExample 3:")
analyze_pattern("5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 0 0 0 0 0 0 0",
                "5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 0 0 0 0 0 0 0 5")

# Analyze test input
print("\nTest Input:")
analyze_pattern("4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4", "")