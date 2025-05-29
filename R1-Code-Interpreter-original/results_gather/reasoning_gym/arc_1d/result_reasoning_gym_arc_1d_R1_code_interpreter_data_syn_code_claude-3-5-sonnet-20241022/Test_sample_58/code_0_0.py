def analyze_pattern(examples):
    for i, (input_str, output_str) in enumerate(examples):
        input_nums = [int(x) for x in input_str.split()]
        output_nums = [int(x) for x in output_str.split()]
        
        # Count zeros
        input_zeros = input_nums.count(0)
        output_zeros = output_nums.count(0)
        
        # Find positions of zeros
        input_zero_pos = [i for i, x in enumerate(input_nums) if x == 0]
        output_zero_pos = [i for i, x in enumerate(output_nums) if x == 0]
        
        print(f"\nExample {i+1}:")
        print(f"Input zeros: {input_zeros} at positions {input_zero_pos}")
        print(f"Output zeros: {output_zeros} at positions {output_zero_pos}")
        print(f"Zero shift: {[out-inp for inp, out in zip(input_zero_pos, output_zero_pos)]}")

examples = [
    ("5 5 0 0 0 0 5 5 5 5 5 5", "0 0 0 5 5 5 5 5 5 5 5 0"),
    ("8 8 8 8 8 8 8 0 8 8 8 8", "8 8 8 8 0 8 8 8 8 8 8 8"),
    ("5 5 5 5 5 5 5 5 0 0 5 5", "5 5 5 5 5 0 0 5 5 5 5 5")
]

analyze_pattern(examples)

# Now analyze the test input
test_input = "6 6 6 6 6 6 6 6 0 0 0 0"
print("\nTest input analysis:")
test_nums = [int(x) for x in test_input.split()]
test_zeros = test_nums.count(0)
test_zero_pos = [i for i, x in enumerate(test_nums) if x == 0]
print(f"Zeros: {test_zeros} at positions {test_zero_pos}")