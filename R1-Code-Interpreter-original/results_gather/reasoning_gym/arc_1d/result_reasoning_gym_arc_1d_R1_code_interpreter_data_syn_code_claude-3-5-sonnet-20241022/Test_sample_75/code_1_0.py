def find_pattern(input_grid):
    # Convert input to list
    nums = [int(x) for x in input_grid.split()]
    
    # Find the start of the main sequence
    # Look for longest stable sequence in middle
    best_start = 0
    max_len = 0
    
    for i in range(len(nums)):
        # Check sequence starting at i
        non_zero_count = sum(1 for j in range(i, len(nums)) if nums[j] != 0)
        if non_zero_count > max_len:
            max_len = non_zero_count
            best_start = i
            
        # If we found a good sequence with multiple non-zero numbers
        if non_zero_count > 10:  # threshold for a "good" sequence
            break
    
    # Create output by moving first part to end
    result = nums[best_start:] + nums[:best_start]
    
    return ' '.join(map(str, result))

# Test input
test_input = "5 5 9 7 4 1 5 6 9 2 7 8 1 7 1 6 8 9 9 8 8 5 4 7 4 9 6 7 7"
print(find_pattern(test_input))