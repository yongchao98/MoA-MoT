def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = [0] * len(input_nums)
    n = len(input_nums)
    
    # First identify sequences of 3 or more
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Check if this is start of a sequence of 3 or more
            count = 1
            j = i + 1
            while j < n and input_nums[j] == input_nums[i]:
                count += 1
                j += 1
            
            if count >= 3:
                # Copy sequence as is
                for k in range(i, i + count):
                    output[k] = input_nums[k]
                i += count
            else:
                # Single number: copy it and add exactly 3 copies immediately after
                num = input_nums[i]
                output[i] = num
                for k in range(1, 4):
                    if i + k < n and input_nums[i + k] == 0:
                        output[i + k] = num
                i += 1
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))