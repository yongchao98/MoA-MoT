def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = [0] * len(input_nums)
    n = len(input_nums)
    
    # First identify and copy sequences of length 3 or more
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Check if this is start of a sequence
            seq_len = 1
            j = i + 1
            while j < n and input_nums[j] == input_nums[i]:
                seq_len += 1
                j += 1
            
            if seq_len >= 3:
                # Copy the sequence as is
                for k in range(i, i + seq_len):
                    output[k] = input_nums[k]
                i += seq_len
            else:
                # Single number: copy it and 3 more times if possible
                num = input_nums[i]
                output[i] = num  # Original number
                for k in range(i + 1, min(i + 4, n)):
                    if k < n and (input_nums[k] == 0 or k > j):
                        output[k] = num
                i += 1
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))