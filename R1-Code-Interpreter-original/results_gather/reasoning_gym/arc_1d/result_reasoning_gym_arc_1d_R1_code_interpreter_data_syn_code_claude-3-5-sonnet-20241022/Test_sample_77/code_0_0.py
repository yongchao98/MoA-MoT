def verify_pattern():
    # Example pairs
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
        
        # Verify if output is input shifted left with 0 at end
        expected_output = input_list[1:] + ['0']
        actual_output = output_list
        
        print(f"Example {i+1}:")
        print(f"Expected: {' '.join(expected_output)}")
        print(f"Actual  : {' '.join(actual_output)}")
        print(f"Match: {expected_output == actual_output}\n")

    # Apply rule to test input
    test_input = "0 3 7 9 4 0 5 0 0 0 0 2 0 3 2 3 0"
    test_input_list = test_input.split()
    predicted_output = ' '.join(test_input_list[1:] + ['0'])
    print(f"Predicted output for test input:")
    print(predicted_output)

verify_pattern()