def find_solution(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    
    # Find if there's a sequence at start
    start_seq_end = 0
    for i in range(n):
        if grid[i] != 0:
            start_seq_end = i + 1
        elif start_seq_end > 0:
            break
    
    if start_seq_end > 0:
        # There's a sequence at start
        start_seq = grid[:start_seq_end]
        if len(start_seq) > 4:
            # Move last 4 to start, rest to end
            output = start_seq[-4:] + [0] * (n - len(start_seq)) + start_seq[:-4]
        else:
            # Keep sequence at start
            output = start_seq + [0] * (n - len(start_seq))
    else:
        # No sequence at start, shift everything left by 3
        output = grid[-3:] + grid[:-3]
    
    return " ".join(map(str, output))

# Test with the input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print(find_solution(test_input))