def process_grid(input_grid):
    input_list = input_grid.split()
    input_nums = [int(x) for x in input_list]
    output = input_nums.copy()  # Start with a copy of input
    n = len(input_nums)
    
    # First pass: Mark positions that are part of sequences (3 or more)
    in_sequence = [False] * n
    i = 0
    while i < n:
        if input_nums[i] != 0:
            # Count consecutive same numbers
            count = 1
            j = i + 1
            while j < n and input_nums[j] == input_nums[i]:
                count += 1
                j += 1
            
            # If sequence of 3 or more found, mark all positions
            if count >= 3:
                for k in range(i, i + count):
                    in_sequence[k] = True
            i = j
        else:
            i += 1
    
    # Second pass: Process isolated numbers
    for i in range(n):
        if input_nums[i] != 0 and not in_sequence[i]:
            num = input_nums[i]
            copies_needed = 3
            j = i + 1
            
            # Make exactly 3 copies, skipping over positions in sequences
            while j < n and copies_needed > 0:
                if not in_sequence[j]:
                    output[j] = num
                    copies_needed -= 1
                j += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 7 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print(process_grid(test_input))