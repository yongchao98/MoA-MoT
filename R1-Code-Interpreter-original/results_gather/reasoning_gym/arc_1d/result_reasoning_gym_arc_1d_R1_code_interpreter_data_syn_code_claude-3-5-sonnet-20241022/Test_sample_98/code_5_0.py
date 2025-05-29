def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = input_nums.copy()
    n = len(input_nums)
    
    # First mark positions that are part of sequences
    in_sequence = [False] * n
    for i in range(n):
        if i > 0 and input_nums[i] == input_nums[i-1] and input_nums[i] != 0:
            in_sequence[i] = True
            in_sequence[i-1] = True
        if i < n-1 and input_nums[i] == input_nums[i+1] and input_nums[i] != 0:
            in_sequence[i] = True
            in_sequence[i+1] = True
    
    # Process isolated numbers
    for i in range(n):
        if input_nums[i] != 0 and not in_sequence[i]:
            num = input_nums[i]
            # Make 3 copies after the number, but stop at any non-zero number
            copies_made = 0
            j = i + 1
            while j < n and copies_made < 3 and input_nums[j] == 0:
                output[j] = num
                copies_made += 1
                j += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))