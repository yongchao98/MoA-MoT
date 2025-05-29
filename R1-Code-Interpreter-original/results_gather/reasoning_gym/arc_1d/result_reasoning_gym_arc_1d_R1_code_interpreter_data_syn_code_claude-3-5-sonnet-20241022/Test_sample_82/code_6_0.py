def find_solution(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    output = [0] * n
    
    # Find if there's a sequence at start
    start_seq_end = 0
    for i in range(n):
        if grid[i] != 0:
            start_seq_end = i + 1
        elif start_seq_end > 0:
            break
    
    # Find other sequences
    other_seqs = []
    start = None
    for i in range(start_seq_end, n):
        if grid[i] != 0:
            if start is None:
                start = i
        elif start is not None:
            other_seqs.append((start, i))
            start = None
    if start is not None:
        other_seqs.append((start, n))
    
    if start_seq_end > 0:
        start_seq = grid[:start_seq_end]
        if len(start_seq) > 4:
            # Move last 4 to start, rest to end
            output[:4] = start_seq[-4:]
            output[-len(start_seq[:-4]):] = start_seq[:-4]
        else:
            # Move sequence 2 positions right
            for i, num in enumerate(start_seq):
                output[(i + 2) % n] = num
    
    # Process other sequences
    for start, end in other_seqs:
        seq = grid[start:end]
        new_start = (start - 3) % n
        for i, num in enumerate(seq):
            output[(new_start + i) % n] = num
    
    return " ".join(map(str, output))

# Test with the input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print(find_solution(test_input))