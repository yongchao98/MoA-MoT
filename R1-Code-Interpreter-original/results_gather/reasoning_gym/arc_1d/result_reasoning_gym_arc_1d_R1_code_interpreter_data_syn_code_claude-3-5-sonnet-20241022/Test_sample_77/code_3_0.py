def analyze_zero_sequences():
    examples = [
        ("9 4 0 0 0 0 0 5 3 0 9 0 9 1 0 0 3",
         "9 4 0 0 0 0 5 3 0 9 0 9 1 0 0 3 0"),
        ("0 0 0 7 3 1 7 0 5 0 0 0 0 0 0 8 8",
         "0 0 7 3 1 7 0 5 0 0 0 0 0 0 8 8 0"),
        ("0 0 8 4 7 5 8 0 4 0 0 0 2 4 6 0 0",
         "0 8 4 7 5 8 0 4 0 0 0 2 4 6 0 0 0")
    ]
    
    for i, (input_str, output_str) in enumerate(examples):
        input_list = input_str.split()
        output_list = output_str.split()
        
        print(f"\nExample {i+1}:")
        # Find sequences of zeros in input
        zero_sequences = []
        current_seq = []
        for j, num in enumerate(input_list):
            if num == '0':
                current_seq.append(j)
            elif current_seq:
                zero_sequences.append(current_seq)
                current_seq = []
        if current_seq:
            zero_sequences.append(current_seq)
            
        print(f"Input zero sequences at positions: {zero_sequences}")
        print(f"Input:  {' '.join(input_list)}")
        print(f"Output: {' '.join(output_list)}")

    # Analyze test input
    test_input = "0 3 7 9 4 0 5 0 0 0 0 2 0 3 2 3 0"
    test_list = test_input.split()
    zero_sequences = []
    current_seq = []
    for j, num in enumerate(test_list):
        if num == '0':
            current_seq.append(j)
        elif current_seq:
            zero_sequences.append(current_seq)
            current_seq = []
    if current_seq:
        zero_sequences.append(current_seq)
    
    print("\nTest input zero sequences at positions:", zero_sequences)
    print("Test input:", test_input)

analyze_zero_sequences()