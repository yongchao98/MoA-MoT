def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    def get_nonzero_sequences(arr):
        sequences = []
        current_seq = []
        for i, num in enumerate(arr):
            if num != 0:
                current_seq.append((num, i))
            elif current_seq:
                sequences.append(current_seq)
                current_seq = []
        if current_seq:
            sequences.append(current_seq)
        return sequences

    input_seqs = get_nonzero_sequences(input_list)
    output_seqs = get_nonzero_sequences(output_list)
    
    print(f"Input sequences: {input_seqs}")
    print(f"Output sequences: {output_seqs}")

# Analyze each example
examples = [
    ("0 0 0 0 0 0 0 5 6 3 9 0", "3 9 0 0 0 0 0 0 0 0 5 6"),
    ("0 0 1 5 3 0 0 0 0 0 0 0", "0 0 0 0 0 1 5 3 0 0 0 0"),
    ("8 8 0 0 0 0 0 0 0 0 0 0", "0 0 0 8 8 0 0 0 0 0 0 0")
]

for i, (inp, out) in enumerate(examples, 1):
    print(f"\nExample {i}:")
    analyze_pattern(inp, out)

# Analyze test input
test_input = "6 3 4 2 6 9 8 1 7 0 0 0"
print("\nTest input analysis:")
analyze_pattern(test_input, "0 0 0 0 0 0 0 0 0 0 0 0")  # dummy output for analysis