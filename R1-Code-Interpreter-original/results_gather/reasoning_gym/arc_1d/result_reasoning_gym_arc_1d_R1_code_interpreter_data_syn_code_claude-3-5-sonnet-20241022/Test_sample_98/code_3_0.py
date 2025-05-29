def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = [0] * len(input_nums)
    n = len(input_nums)
    
    # First copy all numbers from input
    for i in range(n):
        output[i] = input_nums[i]
    
    # Process isolated numbers
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Check if it's isolated (not part of a sequence)
            is_isolated = True
            if i > 0 and input_nums[i] == input_nums[i-1]:
                is_isolated = False
            if i < n-1 and input_nums[i] == input_nums[i+1]:
                is_isolated = False
            
            if is_isolated:
                num = input_nums[i]
                # Fill next three positions with the same number
                for j in range(i+1, min(i+4, n)):
                    output[j] = num
                i += 3  # Skip the positions we just filled
            
        i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))