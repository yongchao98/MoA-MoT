def verify_answer(input_str):
    # Convert input to list of integers
    nums = [int(x) for x in input_str.split()]
    
    # Get counts
    count_2 = nums.count(2)  # 4
    count_0 = nums.count(0)  # 16
    count_1 = nums.count(1)  # 3
    
    # Generate and verify output has correct length
    output = [2]*count_2 + [0]*count_0 + [1]*count_1
    print(" ".join(map(str, output)))
    print(f"Length check: {len(output) == len(nums)}")  # Should be 23

test_input = "2 2 1 0 0 0 0 2 0 0 0 0 0 0 1 0 1 0 2 0 0 0 0"
verify_answer(test_input)