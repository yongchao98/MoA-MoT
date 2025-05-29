def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Example inputs and outputs
    examples = [
        ("7 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 7 7", "7 7 7 7 7 7 7 7 7 7 0 0 0 0 7 7 7 7 7 7"),
        ("3 0 0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3", "0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 0 0 0"),
        ("9 9 9 9 9 9 9 0 0 0 0 0 0 0 0 0 0 0 0 0", "9 9 9 0 0 0 0 0 0 0 0 0 0 0 0 0 9 9 9 9")
    ]
    
    # Analyze test input
    test_input = "0 0 0 0 0 0 0 8 8 8 8 8 8 8 8 0 0 0 0 0"
    test_nums = [int(x) for x in test_input.split()]
    
    # Count numbers in test input
    zeros = test_nums.count(0)
    eights = test_nums.count(8)
    
    # Based on pattern analysis, construct output
    # Pattern seems to be: move 3 of second number to each end
    result = [0] * zeros + [8] * (eights-6) + [0] * 0 + [8] * 3 + [8] * 3
    
    print(" ".join(map(str, result)))

# Test with the given input
analyze_pattern("", "")