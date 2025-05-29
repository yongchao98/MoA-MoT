def find_solution(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    output = [0] * n
    
    # Find the rightmost non-zero sequence
    right_seq_start = None
    right_seq_end = None
    
    for i in range(len(grid)-1, -1, -1):
        if grid[i] != 0:
            if right_seq_end is None:
                right_seq_end = i + 1
            right_seq_start = i
        elif right_seq_end is not None and grid[i] == 0:
            break
    
    if right_seq_start is None:
        return " ".join(map(str, output))
    
    # Get the sequence
    sequence = grid[right_seq_start:right_seq_end]
    
    # Shift left by 3
    new_start = (right_seq_start - 3) % n
    
    # Place the sequence in output
    for i, num in enumerate(sequence):
        output[(new_start + i) % n] = num
    
    # If there are any numbers before right_seq_start, wrap them to the end
    if right_seq_start > 0:
        remaining = [x for x in grid[:right_seq_start] if x != 0]
        for i, num in enumerate(remaining):
            output[-(len(remaining)-i)] = num
    
    return " ".join(map(str, output))

# Test with the input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print(find_solution(test_input))