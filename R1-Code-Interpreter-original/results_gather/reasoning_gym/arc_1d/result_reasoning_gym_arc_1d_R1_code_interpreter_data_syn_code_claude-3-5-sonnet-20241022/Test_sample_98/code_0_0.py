def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = input_nums.copy()
    n = len(input_nums)
    
    # Process each position
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Check if it's part of a sequence
            if i > 0 and input_nums[i] == input_nums[i-1]:
                i += 1
                continue
            if i < n-1 and input_nums[i] == input_nums[i+1]:
                i += 1
                continue
                
            # If isolated number, create pattern
            num = input_nums[i]
            # Fill pattern to the right
            j = i
            while j < n:
                if j == i:
                    output[j] = num
                elif (j - i) % 2 == 0:
                    output[j] = num
                else:
                    output[j] = 0
                j += 1
        i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))