def analyze_transformations():
    examples = [
        ("9 4 0 0 0 0 0 5 3 0 9 0 9 1 0 0 3",
         "9 4 0 0 0 0 5 3 0 9 0 9 1 0 0 3 0"),
        ("0 0 0 7 3 1 7 0 5 0 0 0 0 0 0 8 8",
         "0 0 7 3 1 7 0 5 0 0 0 0 0 0 8 8 0"),
        ("0 0 8 4 7 5 8 0 4 0 0 0 2 4 6 0 0",
         "0 8 4 7 5 8 0 4 0 0 0 2 4 6 0 0 0")
    ]
    
    for i, (input_str, output_str) in enumerate(examples):
        input_list = list(map(int, input_str.split()))
        output_list = list(map(int, output_str.split()))
        
        print(f"\nExample {i+1}:")
        print("Input:", input_list)
        print("Output:", output_list)
        
        # Find the first position where they differ
        first_diff = None
        for j in range(len(input_list)):
            if input_list[j] != output_list[j]:
                first_diff = j
                break
        
        if first_diff is not None:
            print(f"First difference at position {first_diff}")
            print(f"Input before diff: {input_list[:first_diff]}")
            print(f"Output before diff: {output_list[:first_diff]}")
    
    # Test the pattern on new input
    test_input = "0 3 7 9 4 0 5 0 0 0 0 2 0 3 2 3 0"
    test_list = list(map(int, test_input.split()))
    print("\nTest input:", test_list)

analyze_transformations()