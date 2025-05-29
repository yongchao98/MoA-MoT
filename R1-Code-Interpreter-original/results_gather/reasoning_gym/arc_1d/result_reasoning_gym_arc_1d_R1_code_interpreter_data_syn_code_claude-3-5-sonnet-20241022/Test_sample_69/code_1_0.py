def apply_rule(input_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [0] * len(input_list)
    
    # Find sequences
    sequences = []
    current_seq = []
    for i, num in enumerate(input_list):
        if num != 0:
            current_seq.append(i)
        elif current_seq:
            sequences.append(current_seq)
            current_seq = []
    if current_seq:
        sequences.append(current_seq)
    
    # Apply rules to each sequence
    for seq_idx, seq in enumerate(sequences):
        value = input_list[seq[0]]
        
        if seq_idx == 0:  # First sequence
            # Copy all but last number
            for i in seq[:-1]:
                output_list[i] = value
        elif seq_idx == len(sequences) - 1:  # Last sequence
            # Copy all numbers and add one at the end
            for i in seq:
                output_list[i] = value
            output_list[-1] = value
        else:  # Middle sequences
            # Copy sequence as is
            for i in seq:
                output_list[i] = value
    
    return ' '.join(map(str, output_list))

test_input = "4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4"
result = apply_rule(test_input)
print(result)