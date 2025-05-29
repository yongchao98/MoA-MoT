def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = [0] * len(input_nums)
    n = len(input_nums)
    
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Check if it's part of a sequence
            is_sequence = False
            if (i > 0 and input_nums[i] == input_nums[i-1]) or \
               (i < n-1 and input_nums[i] == input_nums[i+1]):
                output[i] = input_nums[i]
                is_sequence = True
            
            if not is_sequence:
                num = input_nums[i]
                # If at start or after zeros, create pattern left
                if i == 0:
                    output[0] = num
                    for j in range(1, min(4, n)):
                        output[j] = num
                
                # Create pattern right
                for j in range(i, min(i + 4, n)):
                    output[j] = num
                if i + 4 < n:
                    output[i + 4] = 0
        i += 1
    
    # Preserve original sequences
    for i in range(n):
        if i > 0 and input_nums[i] == input_nums[i-1]:
            output[i] = input_nums[i]
            output[i-1] = input_nums[i-1]
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))