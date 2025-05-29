def analyze_pattern():
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
        # Find where input and output differ
        for j in range(len(input_list)):
            if input_list[j] != output_list[j]:
                print(f"First difference at position {j}")
                print(f"Input subsequence: {' '.join(input_list[j-2:j+3])}")
                print(f"Output subsequence: {' '.join(output_list[j-2:j+3])}")

    # Apply to test input
    test_input = "0 3 7 9 4 0 5 0 0 0 0 2 0 3 2 3 0"
    test_list = test_input.split()
    # Find first non-zero sequence
    zeros = []
    for i, num in enumerate(test_list):
        if num == "0":
            zeros.append(i)
    print("\nTest input zero positions:", zeros)
    print("Test input:", test_input)

analyze_pattern()