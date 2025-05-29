def find_pattern(input_str):
    nums = [int(x) for x in input_str.split()]
    n = len(nums)
    non_zeros = [x for x in nums if x != 0]
    zeros = [x for x in nums if x == 0]
    
    # Find the ideal center position for zeros
    zero_count = len(zeros)
    start_pos = (n - zero_count) // 2
    
    # Construct output
    result = non_zeros[:start_pos] + zeros + non_zeros[start_pos:]
    print(f"Constructed output: {result}")
    return result

# Test with all examples
examples = [
    "5 5 0 0 0 0 5 5 5 5 5 5",
    "8 8 8 8 8 8 8 0 8 8 8 8",
    "5 5 5 5 5 5 5 5 0 0 5 5"
]

for ex in examples:
    print("\nTesting example:", ex)
    find_pattern(ex)

# Test with the input
test = "6 6 6 6 6 6 6 6 0 0 0 0"
print("\nTesting input:", test)
find_pattern(test)