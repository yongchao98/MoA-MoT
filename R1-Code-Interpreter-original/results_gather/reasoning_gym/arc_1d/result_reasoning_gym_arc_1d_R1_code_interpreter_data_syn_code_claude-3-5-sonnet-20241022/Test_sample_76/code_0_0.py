def analyze_pattern(input_str, output_str):
    input_list = input_str.split()
    output_list = output_str.split()
    
    # Count consecutive zeros from start and end
    def count_zeros(lst, from_start=True):
        count = 0
        idx = 0 if from_start else len(lst)-1
        step = 1 if from_start else -1
        while 0 <= idx < len(lst) and lst[idx] == '0':
            count += 1
            idx += step
        return count
    
    # Print analysis
    print(f"Input: {input_list}")
    print(f"Output: {output_list}")
    print(f"Leading zeros: {count_zeros(input_list)} -> {count_zeros(output_list)}")
    print(f"Trailing zeros: {count_zeros(input_list, False)} -> {count_zeros(output_list, False)}")
    print("---")

# Analyze examples
examples = [
    ("5 5 5 5 5 0 5 5 5 5 5 5", "5 5 5 5 5 5 0 5 5 5 5 5"),
    ("0 0 0 2 2 2 2 2 2 2 0 0", "0 0 0 0 2 2 2 2 2 2 2 0"),
    ("0 0 0 3 3 3 3 3 0 0 0 0", "0 0 0 0 3 3 3 3 3 0 0 0")
]

for input_str, output_str in examples:
    analyze_pattern(input_str, output_str)

# Analyze test case
test_input = "6 6 6 0 0 0 0 0 0 0 0 6"
print("Test case analysis:")
analyze_pattern(test_input, "")