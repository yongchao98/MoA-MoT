def find_longest_sequence(input_str):
    numbers = [int(x) for x in input_str.split()]
    output = [0] * len(numbers)
    
    # Find all sequences
    current_seq = []
    longest_seq = []
    
    for i, num in enumerate(numbers):
        if num != 0:
            current_seq.append(i)
        elif current_seq:
            if len(current_seq) > len(longest_seq):
                longest_seq = current_seq.copy()
            current_seq = []
            
    # Check last sequence
    if current_seq and len(current_seq) > len(longest_seq):
        longest_seq = current_seq
    
    # Keep only longest sequence
    for i in longest_seq:
        output[i] = numbers[i]
    
    return ' '.join(map(str, output))

test_input = "0 5 0 5 0 0 0 0 0 0 5 0 0 5 5 5 0 0 0 0"
result = find_longest_sequence(test_input)
print(result)