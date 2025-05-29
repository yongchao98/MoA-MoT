def analyze_continuous_sequences(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find continuous sequences of non-zero numbers
    def find_sequences(arr):
        sequences = []
        current_seq = []
        for i, num in enumerate(arr):
            if num != 0:
                current_seq.append(i)
            elif current_seq:
                sequences.append(current_seq)
                current_seq = []
        if current_seq:
            sequences.append(current_seq)
        return sequences
    
    in_seq = find_sequences(input_list)
    out_seq = find_sequences(output_list)
    
    print("Input sequences:", in_seq)
    print("Output sequences:", out_seq)
    print("Length of longest input sequence:", max([len(s) for s in in_seq]) if in_seq else 0)
    print("Length of longest output sequence:", max([len(s) for s in out_seq]) if out_seq else 0)

print("Example 1:")
analyze_continuous_sequences("0 0 7 7 7 7 7 7 7 7 7 0 7 0 0 7 0 0 7 0", 
                          "0 0 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0 0 0")

print("\nExample 2:")
analyze_continuous_sequences("1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1", 
                          "0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1")

print("\nExample 3:")
analyze_continuous_sequences("0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0", 
                          "0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0")