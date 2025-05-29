def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = input_nums.copy()  # Start with exact copy
    n = len(input_nums)
    
    # First identify all positions that are part of sequences of 3 or more
    in_sequence = [False] * n
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Count consecutive numbers
            start = i
            count = 1
            while i + count < n and input_nums[i + count] == input_nums[i]:
                count += 1
            
            # Mark if it's a sequence of 3 or more
            if count >= 3:
                for j in range(start, start + count):
                    in_sequence[j] = True
            
            i += count
        else:
            i += 1
    
    # Now process isolated numbers
    i = 0
    while i < n:
        if input_nums[i] != 0 and not in_sequence[i]:
            # Check next three positions
            num = input_nums[i]
            for j in range(1, 4):
                if i + j < n and not in_sequence[i + j]:
                    output[i + j] = num
        i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))